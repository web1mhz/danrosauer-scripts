# Dan Rosauer September 2012

##STEPS
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

import arcpy
import sys, os, math, numpy, csv, os, string
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

sys.path.append("C:\\Users\\u3579238\Work\AMT\\")
import LineageFunctions

### PARAMETERS ###
higher_taxon = "geckoes"
genus_list = ["Heteronotia"]

base_dir = "C:\\Users\\u3579238\\work\\AMT\\Models\\"
target_location = base_dir + "lineage_models\\"  # where the lineage model grids and working data
output_gdb_name = "lineage_models.gdb"
scratch_workspace = "C:\\Users\\u3579238\\Work\\AMT\\Models\\"  #scratch workspace is used by ArcGIS for temporary files during analysis
scratch_gdb_name = "scratch.gdb"
export_asc = True
asc_target_location = base_dir + "lineage_models\\asc_cube_method\\"

buffer_dist = 2.5      # the buffer distance in degrees
additional_buffer = 0  ## how much (as a proportion) the output grids should extend beyond the buffered points
grid_resolution = 0.01
spRef = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
Australia_extent = arcpy.Extent(112.9,-43.75,153.64,-9)
Lineage_field_name = "Lineage"
Distance_method = "model-cost"      ## determines whether distance is calculated as euclidean or model-weighted cost distance
                                    ## so far, can be "euclidian" or "model-cost"
Weight_function = "inverse_cube"    ## determines whether lineage weight is calculated as 1/distance or 1/(distance^2), or simply closest distance
                                    ## so far, can be "inverse" or "inverse_square" or "cost_allocation"
Min_dist_value = 0.0001             ## sets a floor for cost. This prevents a division by 0 error.  In future change minimum value to grid_resolution / 2 (distance to edge of cell)
                                    ##   but keep current value for consistency in this study
Min_weight_threshold = 0.02         ## weights below this for any layer are set to 0.  If the value here is 0, then no threshold is applied
Scale_to = "model"                  ## determines whether lineage weights sum to the model suitability or to 1
                                    ## can be "model" or "one"
model_edited_suffix = "_edited"     ## this is to check if there is a version of the model cliiped to reduce overprediction

handle_minor = ""        ## if true, then lineages with < than the specified proportion of the lineage sum 
omit_minor_threshold = 0.1         ## for that cell, are set to 0

skip_distance_layers = False         ## skip creating the distance layers - they are already done.  THIS OPTION IS ONLY TO SAVE TIME DURING DEBUGGING

# a changeable list to allow for species in the dataset to be skipped
named_species   = []
use_list        = ""  #specify whether to:
                        #do - the named species (use_list="do")
                        #skip - the named species (use_list="skip")
                        #do all the species in the data and ignore the named species list (use_list="" or anything else);
Lin_exclude_list = ["_","lin_","CZ","lin_CZ","0"]
Lin_exclude_list = Lin_exclude_list + ["johnstonei_3"]  # this is a temporary workaround for Carlia johnstonei lineage 3, which crashes the script
# because it is on Kimberley islands beyond the raster environmental data

# HETERONOTIA BINOEI SPECIFIC CHANGES
Australia_extent = arcpy.Extent(112.9,-22,153.64,-9) # modified for H.binoei study
lineage_site_filename = base_dir + "lineage_sites\\" + "Heteronotia_lin_loc_from_db_18mar14_nth_of_22S.csv"
maxent_model = base_dir + "species_models\\maxent\\maxent_models.gdb" + "\\Heteronotia_binoei_17mar14_mean"

for genus in genus_list:
    print "\nGenus: " + genus + "\n"

    #lineage_site_filename = base_dir + "lineage_sites\\" + genus + "_lin_loc.csv"
    #maxent_model_base = base_dir + "species_models\\maxent\\maxent_models.gdb"
        
    # create the scratch geodatabase if needed (to store temporary layers during analysis)
    try:
        arcpy.CreateFileGDB_management(scratch_workspace, scratch_gdb_name)
    except:
        print "File geodatabase "+ scratch_workspace + " " + scratch_gdb_name + " already exists or creation failed"
    
    # create the output geodatabase if needed
    if not os.path.exists(target_location):
        os.makedirs(target_location)
    try:
        target_location_ESRI = string.replace(target_location,"\\","/")
        arcpy.CreateFileGDB_management(target_location_ESRI,output_gdb_name)
    except:
        print "File geodatabase " + target_location + output_gdb_name + " already exists or creation failed"
    
    # Load the sequence site data
    print "\nLoading the sequenced sites\n"
    
    # make each column being used, into a list
    # and also create lists of unique AnalysisGroups and Lineages
    with open(lineage_site_filename, 'rb') as csvfile:
        sequence_csv = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_MINIMAL)
        rownum = 0
        #rowdata = []
        Species = []
        ModelGroup=[]
        Lineage = []
        Lat=[]
        Long=[]
        GroupList=[]
        GroupLineageList =[]
        UseCol = -1
        
        for row in sequence_csv:
            # Save header row.
            if rownum == 0:
                header = row
                rownum += 1
                
                # find the column called 'use'
                columns = range(len(header))
                for k in columns:
                    if str.lower(header[k]) == "use":
                        UseCol = k
                        break

            else:
                if (UseCol == -1 or row[UseCol] == '1'):
                    try:        # if the lat or long can't be converted to a number, then skip that row
                        Lat.append(float(row[5]))
                        try:
                            Long.append(float(row[6]))
                            # code gets to here for valid lat and long, so other steps can go here too
                            #rowdata.append(row)
                            Species.append(row[2])
                            ModelGroup.append(row[3])
                            Lineage.append(row[4])                    
                            if row[3] not in GroupList:
                                GroupList.append(row[3])
                            if row[1:5] not in GroupLineageList and len(row[4]) > 0:
                                GroupLineageList.append(row[1:5])
                            rownum += 1
                        except:
                            pass
                    except:
                        pass
        # so now we have a list for each column, excluding where 'Use' is not 1 and rows with null coordinates
    
    # set the geoprocessing environment
    env.workspace  = target_location_ESRI + output_gdb_name
    env.scratchWorkspace = scratch_workspace + scratch_gdb_name
    env.extent = Australia_extent
    
    # Make XY Event Layer
    try:
        points_layer = "sequenced_sites"
        #arcpy.MakeXYEventLayer_management(lineage_site_filename, "long", "lat", target_location_ESRI+points_layer, spRef)
        layer_result = arcpy.MakeXYEventLayer_management(lineage_site_filename, "long", "lat", points_layer, spRef)
        #print arcpy.GetCount_management(layer_result) + " sites loaded\n"
        points_layer = layer_result[0]
        
    except:
       # If an error occurred print the message to the screen
        print arcpy.GetMessages()
    
    # export the layer as a feature class
    outLocation = env.workspace
    outFeatureClass = genus + "_lin_points"
    
    try:
        arcpy.env.overwriteOutput=True
        arcpy.CopyFeatures_management(points_layer, outFeatureClass, "", "0", "0", "0")
    except:
        print "\n" + outFeatureClass + " could not be created, or already exists."
    points_fc = outLocation + "/" + outFeatureClass
    
    # restrict the GroupList to particular species based on the names_species parameter
    if use_list == "do":
        GroupList = list(set(named_species).intersection(set(GroupList)))
    elif use_list == "skip":
        GroupList = list(set(GroupList).difference(set(named_species)))
    
    # Loop through the list of Model Groups
    for group in GroupList:
        if group != 0 and group != "":
        
            print "\nStarting group " + group + "\n"
            
            # check for an edited verion of the Maxent model, and use if it exists
            #maxent_model = maxent_model_base + "\\" + string.replace(group," ","_")
            if arcpy.Exists(maxent_model + model_edited_suffix):
                maxent_model = maxent_model + model_edited_suffix
        
            # set the environment
            env.snapRaster = maxent_model
            env.mask = maxent_model
            env.extent=maxent_model
        
            # define spatial data for this group
            maxent_raster = arcpy.sa.Raster(maxent_model)
            groupDefQuery = "[ModelGroup] = '" + group + "'"
            points_layer.definitionQuery = groupDefQuery
        
            # start a list of layers to delete at the end
            layers_to_delete = []
    
            ## get a list of the lineages in this group
            lineage_list=[]
            print "Lineages in " + group + ":"
            for row in GroupLineageList:
                if row[3] not in lineage_list:
                    if row[2] == group and "," not in row[3] and row[3] not in Lin_exclude_list:  # exclude particular lineage names and those with dodgy punctuation (fix in data later)
                        
                        lineage = row[3]
                        lineage = lineage.strip()
                        lineage_list.append(lineage)
                        print "   ", row[3]

            maxent_extent = maxent_raster.extent
            
            # get the extent of the points for the model group
            yrange=LineageFunctions.getFieldMinMax(points_fc,"Lat")
            ymin=yrange[0]
            ymax=yrange[1]
            xrange=LineageFunctions.getFieldMinMax(points_fc,"Long")
            xmin=xrange[0]
            xmax=xrange[1]
    
            buffer_ratio = (1 + additional_buffer)
            extent_buffer = buffer_dist * buffer_ratio
        
            # new extent is the same as points layer + a buffer
            # but where the extended buffer goes beyond the extent of the maxent model, limit to the model extent (determined by the buffered environemnt grids)
            xmin=math.floor(max([xmin - extent_buffer,maxent_extent.XMin]))
            ymin=math.floor(max([ymin - extent_buffer,maxent_extent.YMin]))
            xmax=math.ceil(min([xmax + extent_buffer,maxent_extent.XMax]))
            ymax=math.ceil(min([ymax + extent_buffer,maxent_extent.YMax]))
            env.extent = arcpy.Extent(xmin,ymin,xmax,ymax)
        
            if len(lineage_list) > 1:
            
                ### generate a weight grid for each lineage  START OF STEP 4
                
                ## get a selectable layer for the sequenced sites
                lin_lyr = arcpy.MakeFeatureLayer_management(points_fc,"lineage_layer")[0]
            
                print "\nLooping through the lineages in " + group + " to generate weight grids\n"
                count = 0
                
                if Distance_method == "model-cost":
                    model_cost = -1 * arcpy.sa.Ln(maxent_model)
    
                if not skip_distance_layers:
                
                    for lineage in lineage_list:
                        count += 1
                    
                        where_clause = '"ModelGroup" = ' + "'" + group + "' and " + '"' + Lineage_field_name + '"' + " = " + "'" + lineage + "'"
                        arcpy.SelectLayerByAttribute_management(lin_lyr, "NEW_SELECTION", where_clause)
                        arcpy.CopyFeatures_management(lin_lyr, "lineage_points", "", "0", "0", "0")
                        layers_to_delete.append("lineage_points")
                    
                        # create a distance layer for the current lineage
                        if str(lineage) == "0":
                            print "Creating distance layer for sequenced locations of unknown lineage"
                        else:
                            print "Creating distance layer for lineage " + lineage
                        lineage_weight_gridname = "lin_weight_" + string.replace(group," ","_") +"_" + lineage
                        lineage_weight_gridname = lineage_weight_gridname.rstrip()
                        
                        if Distance_method == "model-cost":                                   ## STEP 5b
                            ## calculates the least cost distance to the nearest lineage point
                            ## the result is written directly to lineage_dist_gridname
                            lin_dist = arcpy.sa.PathDistance("lineage_points",model_cost)
                            # change zero values to a very small non-zero value, to avoid nodata in division
                            lin_dist = arcpy.sa.Con(lin_dist==0,0.0001,lin_dist)
                            
                            if Weight_function == "inverse_square":                 ## STEP 6b
                                lin_weight = 1/(lin_dist ** 2)
                            elif Weight_function == "inverse_cube":
                                lin_weight = 1/(lin_dist ** 3)
                            elif Weight_function == "inverse_quad":
                                lin_weight = 1/(lin_dist ** 4)
                            else:
                                lin_weight = 1/lin_dist                             ## STEP 6a       this comes 2nd, as it is the default, for any other values of Weight_function            
                    
                        else:
                            lin_dist = arcpy.sa.EucDistance(lin_lyr,"",grid_resolution)     ## STEP 5a
                            # change zero values to a very small non-zero value, to avoid nodata in division
                            lin_dist = arcpy.sa.Con(lin_dist > Min_dist_value,Min_dist_value,lin_dist)
                            
                            if Weight_function == "inverse_square":                 ## STEP 6b
                                lin_weight = 1/(lin_dist ** 2)
                            elif Weight_function == "inverse_cube":
                                lin_weight = 1/(lin_dist ** 3)
                            elif Weight_function == "inverse_quad":
                                lin_weight = 1/(lin_dist ** 4)
                            else:
                                lin_weight = 1/lin_dist                             ## STEP 6a       this comes 2nd, as it is the default, for any other values of Weight_function            
                       
                        # reset points definition to whole Model Group
                        arcpy.SelectLayerByAttribute_management(lin_lyr, "CLEAR_SELECTION")
                    
                        # apply a threshold to each weight grid                 ## STEP 7
                        if Min_weight_threshold > 0:
                            where_clause = '"VALUE" >= ' + str(Min_weight_threshold)
                            lin_weight = arcpy.sa.Con(lin_weight, lin_weight, 0, where_clause)
                            
                        # set NoData values to 0
                        lin_weight=arcpy.sa.Con(arcpy.sa.IsNull(lin_weight),0,lin_weight)
                    
                        lin_weight.save(lineage_weight_gridname)  ## this layer should be kept until the final weights are calculated
                        
                        #lineage_weight_gridname = arcpy.sa.Con(lineage_weight_gridname, 0, lineage_weight_gridname, lineage_weight_gridname = NULL)
    
                        layers_to_delete.append(lineage_weight_gridname)
                    
                        # calculate the sum of weights for each pixel to scale values later  ## STEP 8
                        if count == 1:
                            weight_sum = lin_weight
                        else:
                            weight_sum = weight_sum + lin_weight
                    
                    print "\nDistance layers done.\n"
                
                else:
                    print "\nDistance layer creation skipped because skip_distance_layers = True.\nThis should only be used during debugging.\n"
                
               # remove very small values for each pixel and rescale so the rest add up to the species model
                if handle_minor == "threshold":
                    count = 0
                    for lineage in lineage_list:                                ## STEPS 9 and 10
                        count += 1
                        if str(lineage) != "0":  #lineage 0 is used to refer to sequenced locations without a named lineage.
                            print "Removing lineage values which account for less than " + (str(omit_minor_threshold*100)) + "% of pixel score. " + lineage
                            lineage_weight_gridname = "lin_weight_" + string.replace(group," ","_") +"_" + lineage 
                            lin_weight = arcpy.sa.Raster(lineage_weight_gridname)
                            
                            lin_proportion = lin_weight / weight_sum
                            where_clause = "VALUE > " + str(omit_minor_threshold)
                            lin_weight = arcpy.sa.Con(lin_proportion,lin_weight,0,where_clause)
                            lin_weight.save(lineage_weight_gridname)
                        if count == 1:
                            new_weight_sum = lin_weight
                        else:
                            new_weight_sum = new_weight_sum + lin_weight
                    weight_sum = new_weight_sum
                    layers_to_delete.append(lin_proportion)
                    layers_to_delete.append(new_weight_sum)
                    print "\nFinished removing low scores\n"

                # calculate the scaled weight for each lineage, so they sum to a) 1 or b) the model suitability
                for lineage in lineage_list:                                ## STEPS 9 and 10
                    if str(lineage) != "0":  #lineage 0 is used to refer to sequenced locations without a named lineage.
                        print "creating final scaled weight grid for lineage " + lineage
                        lineage_weight_gridname = "lin_weight_" + string.replace(group," ","_") +"_" + lineage 
                        lin_weight = arcpy.sa.Raster(lineage_weight_gridname)
                    
                        if Scale_to == "model":
                            lin_weight_scaled = (lin_weight / weight_sum) * maxent_raster
                        else:
                            lin_weight_scaled = (lin_weight / weight_sum)
                    
                        lineage_scaled_weight_name = "lin_model_" + string.replace(group," ","_") +"_" + lineage
                        lin_weight_scaled.save(lineage_scaled_weight_name)  ## THIS is a final layer to keep!
                        
                        if export_asc:
                            asc_filename = asc_target_location + lineage_scaled_weight_name + ".asc"
                            arcpy.RasterToASCII_conversion(lin_weight_scaled, asc_filename)
                            arcpy.DefineProjection_management(asc_filename, spRef)
                
                    else:   # lineage 0 is not an actual lineage, but is used to refer to sequenced locations without a named lineage.
                            # a final model is not created for lineage 0, but it should affect the values for the other lineages.
                        print "Sequenced locations of unnamed lineage were not modelled, but do affect values for other lineages."
            
            # where the model group has just a single lineage, simply copy the maxent model as the lineage model
            else:
                lineage = lineage_list[0]
                lineage_scaled_weight_name = "lin_model_" + string.replace(group," ","_") +"_" + lineage
                maxent_raster.save(lineage_scaled_weight_name)  # saving the maxent model for the model_group as the lineage model
                print "created lineage grid for lineage " + lineage + " as a copy of the model for group: " + group
                        
                if export_asc:
                    asc_filename = asc_target_location + lineage_scaled_weight_name + ".asc"
                    arcpy.RasterToASCII_conversion(maxent_raster, asc_filename)
                    arcpy.DefineProjection_management(asc_filename, spRef)                
                
            print "\nAnalysis for " + group + " completed - now deleting temporary data."

            # and finally, delete temporary layers
            #layers_to_delete = list(set(layers_to_delete))  # remove duplicates
            for layer in layers_to_delete:
                try:
                    arcpy.Delete_management(layer)
                    print layer + " deleted"
                except:
                    try:
                        print layer + " NOT deleted"
                    except: print " NOT deleted"
        
        print "\n   ****************\n   * FINISHED Model Group :",group ,"*\n   ****************"
                
    env.workspace = scratch_workspace + scratch_gdb_name
    datasets = arcpy.ListDatasets()
    for dataset in datasets:
        try:
            arcpy.Delete_management(dataset)
            print dataset + " deleted"
        except:
            print dataset + " NOT deleted"
            
    print "\nLineage models done for " + genus + "\n"
        
print "\nLineage models done!\n"
