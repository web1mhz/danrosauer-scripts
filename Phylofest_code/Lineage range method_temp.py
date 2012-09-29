# Dan Rosauer September 2012

##STEPS TO ENCODE
## 1. import the points for the whole species
##
## 2. to bound the whole analysis, use euc distance to create a grid to define a boundary at a set distance
##
## 3. outside of this code (or from within if we can call Maxent like from R) generate a species distribution model for the whole species
##
## 4. loop through all of the lineages in the species
##    5a. generate a euc distance layer from sequenced locations for each lineage, bounded by the total species euc distance layer from (2)
##     or 
##    5b. generate a cost distance layer from sequenced locations for each lineage, using the maxent model to define the cost.  Cost = 1 - suitability
##    
##    6a. generate a weight layer for each lineage as 1 / distance  from (5)
##     or
##    6b. generate a weight layer for each lineage as 1 / distance^2  from (5)
##
##    7.  set all weights below a threshold to 0, to reduce the effect of distant lineages
##
## 8. sum all of the lineage weight layers
##
## 9. divide each lineage weight layer by the sum of weights (8) so that the weights for each pixel sum to 1
##
## 10. multiply each lineage weight layer by the model likelihood so that the weights for each pixel sum to the model likelihood.

## All of the above could be nested within a loop that iterates through species.

import arcpy, sys, os
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

sys.path.append("C:\\Users\\u3579238\Work\Phylofest\\")
import LineageFunctions

# PARAMETERS ###
species_name = "Diporiphora_bilineata_group"   # of course this would not be a poarameter in the multi-species version
points = "C:\\Users\\u3579238\Work\Phylofest\\Diporiphora_seq.shp"  # sequenced locations point shapefile
output_location = "C:\\Users\\u3579238\Work\Phylofest\\Diporiphora_group.gdb"
maxent_model = "C:\\Users\\u3579238\Work\Phylofest\\Diporiphora_group.gdb\\D_bilineata_grp_maxent"
buffer_dist = 5                                                               # the buffer distance in degrees
additional_buffer = 0  ## how much (as a proportion) the output grids should extend beyond the buffered points
grid_resolution = 0.01
Australia_extent = arcpy.Extent(112,-44,160,-8)
Lineage_field_name = "Lineage"
Distance_method = "euclidian"       ## determines whether distance is calculated as euclidean or model-weighted cost distance
                                    ## so far, can be "euclidian" or "model-cost"
Weight_function = "inverse_square"  ## determines whether lineage weight is calculated as 1/distance or 1/(distance^2)
                                    ## so far, can be "inverse" or "inverse_square"
Min_weight_threshold = 0            ## weights below this for any layer are set to 0.  If the value here is 0, then no threshold is applied
Scale_to = "model"                  ## determines whether lineage weights sum to the model suitability or to 1

print "\nStarting " + species_name + "\n"

                                    ## can be "model" or "one"
# set the environment
env.workspace  = output_location
env.snapRaster = maxent_model
env.extent = Australia_extent

maxent_model = arcpy.sa.Raster(maxent_model)

# start a list of layers to delete at the end
layers_to_delete = []

## Buffer the full set of sequenced locations   STEPS 1 & 2
print "Buffering all points\n"
arcpy.env.mask = "maxent_model"
point_buffer = arcpy.sa.EucDistance(points,buffer_dist,grid_resolution)
point_buffer.save("all_points_buf")  ## once the code is working, no need to save this layer
layers_to_delete.append("all_points_buf")

## from the points shapefile, get a list of the lineages
print "Lineage list\n"
lineage_list = LineageFunctions.getFieldValues(points,Lineage_field_name)
print lineage_list

#calculate the extent for weight grids
points_properties = arcpy.Describe(points)
points_extent = points_properties.extent
extent_buffer = (1 + additional_buffer)
new_extent = [points_extent.xmin - extent_buffer, points_extent.ymin - extent_buffer, points_extent.xmax + extent_buffer, points_extent.ymax + extent_buffer]
env.extent = arcpy.Extent(new_extent) # set the analysis extent

## generate a weight grid for each lineage  START OF STEP 4
arcpy.env.mask = "all_points_buf"
arcpy.MakeFeatureLayer_management(points,"lin_lyr")
layers_to_delete.append("lin_lyr")
print "Looping through the lineages in " + species_name + " to generate weight grids\n"
count = 0

for lineage in lineage_list:
    count = count + 1

    # select points for the current lineage
    where_clause = '"' + Lineage_field_name + '" = \'' + lineage + '\''
    arcpy.SelectLayerByAttribute_management("lin_lyr", "NEW_SELECTION", where_clause)

    # create a distance layer for the current lineage
    print "creating distance layer for lineage " + lineage
    lineage_dist_gridname = "dist_lineage_" + lineage
    
    if Distance_method == "model-cost":                                   ## STEP 5b
        ## calculates the least cost distance to the nearest lineage point
        ## the result is written directly to lineage_dist_gridname
        temp = arcpy.sa.PathAllocation(in_source_data="lin_lyr", in_cost_raster=maxent_model, out_distance_raster=lineage_dist_gridname)    
        lin_dist = arcpy.sa.Raster(lineage_dist_gridname)
        layers_to_delete.append(lineage_dist_gridname)

    else:
        lin_dist = arcpy.sa.EucDistance("lin_lyr","",grid_resolution)     ## STEP 5a
        lin_dist.save(lineage_dist_gridname)  ## once the code is working, no need to save this layer
        layers_to_delete.append(lineage_dist_gridname)

    # unselect lineage points
    arcpy.SelectLayerByAttribute_management("lin_lyr", "CLEAR_SELECTION")

    # create a weight layer for the current lineage
    print "creating weight layer for   lineage " + lineage    
    lineage_weight_gridname = "weight_lineage_" + lineage
    if Weight_function == "inverse_square":                 ## STEP 6b
        lin_weight = 1/(lin_dist ** 2)
    else:
        lin_weight = 1/lin_dist                             ## STEP 6a       this comes 2nd, as it is the default, for any other values of Weight_function

    # apply a threshold to each weight grid                 ## STEP 7
    if Min_weight_threshold > 0:
        lin_weight = Con(lin_weight, lin_weight, 0, "VALUE >= Min_weight_threshold")

    lin_weight.save(lineage_weight_gridname)  ## this layer should be kept until the final weights are calculated
    layers_to_delete.append(lineage_weight_gridname)

    # calculate the sum of weights for each pixel to scale values later  ## STEP 8
    if count == 1:
        weight_sum = lin_weight
    else:
        weight_sum = weight_sum + lin_weight

# calculate the scaled weight for each lineage, so they sum to a) 1 or b) the model suitability
for lineage in lineage_list:                                ## STEPS 9 and 10
    print "creating final scaled weight grid for lineage " + lineage
    lineage_weight_gridname = "weight_lineage_" + lineage
    lin_weight = arcpy.sa.Raster(lineage_weight_gridname)

    if Scale_to == "model":
        lin_weight_scaled = (lin_weight / weight_sum) * maxent_model
    else:
        lin_weight_scaled = (lin_weight / weight_sum)

    lineage_scaled_weight_name = "scaled_weight_lin" + lineage
    lin_weight_scaled.save(lineage_scaled_weight_name)  ## THIS is a final layer to keep!

# and finally, delete temporary layers
print "\nAnalysis completed - now deleting temporary data."
for layer in layers_to_delete:
    arcpy.Delete_management(layer)
    print layer + " deleted"
                    
