# Dan Rosauer November 2012

# a quickly edited version of the lineage range method, to read in the species and lineage distributions and export a .kml for each species.

import arcpy, sys, os, math, numpy, csv, os, string
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

sys.path.append("C:\\Users\\u3579238\Work\Phylofest\\")
import LineageFunctions

### PARAMETERS ###
genus = "Eulamprus"  # genus could refer to any group being handled as a set
higher_taxon = "skinks"
base_dir = "C:\\Users\\u3579238\\work\\Phylofest\\Models\\" + higher_taxon + "\\"
sequence_site_filename = base_dir + "sequence_sites\\" + genus + "_lin_loc.csv"
species_site_filename = base_dir + "species_sites\\" + genus + "_maxent.csv"
target_location = base_dir + "species_sites\\" + "KML\\"  # where the lineage model grids and working data

Australia_extent = arcpy.Extent(112.9,-43.75,153.64,-9)
Lineage_field_name = "Lineage"
 
# Load the sequence site data
print "\nLoading the sequenced sites\n"

# make each column being used, into a list
# and also create lists of unique Model Groups 
with open(sequence_site_filename, 'rb') as csvfile:
    sequence_csv = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_MINIMAL)
    rownum = 0
    GroupList=[]
    Lat=[]
    Long=[]
    
    for row in sequence_csv:
        # Skip header row.
        if rownum == 0:
            rownum += 1
        else:
            try:        # if the lat or long can't be converted to a number, then skip that row
                Lat.append(float(row[5]))
                try:
                    # code gets to here for valid lat and long, so other steps can go here too
                    Long.append(float(row[6]))
                    if row[3] not in GroupList and len(row[3])>0:
                        GroupList.append(row[3])
                except:
                    pass  # a nothing action to meet the syntax rules
            except:
                pass  # a nothing action to meet the syntax rules
    # so now we have a list of model groups, excluding rows with null coordinates

# set the geoprocessing environment
target_location_ESRI = string.replace(target_location,"\\","/")
env.workspace  = target_location_ESRI
env.extent = Australia_extent

# Make XY Event Layer
try:
    spRef = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
    points_layer_lin = "sites"
    #arcpy.MakeXYEventLayer_management(sequence_site_filename, "long", "lat", target_location_ESRI+points_layer_lin, spRef)
    layer_result = arcpy.MakeXYEventLayer_management(sequence_site_filename, "long", "lat", points_layer_lin, spRef)
    points_layer_lin = layer_result[0]
    
except:
   # If an error occurred print the message to the screen
    print arcpy.GetMessages()
    
# export the layer as a feature class
outLocation = env.workspace
outFeatureClass = genus + "_lin_points"

try:
    arcpy.env.overwriteOutput=True
    arcpy.CopyFeatures_management(points_layer_lin, outFeatureClass, "", "0", "0", "0")
except:
    print "\n" + outFeatureClass + " could not be created, or already exists."
points_fc = outLocation + "/" + outFeatureClass

# Loop through the list of Model Groups
for group in GroupList:
    if group != 0 and group != "":
    #if group == "Saproscincus challengeri":   # this is a temp line to model just one group!!
    
        print "\nStarting group " + group + "\n"
        
        groupDefQuery = "[ModelGroup] = '" + group + "'"
        points_layer_lin.definitionQuery = groupDefQuery
               
        # get the extent of the points for the model group
        xrange=LineageFunctions.getFieldMinMax(points_layer_lin,"long")
        xmin=xrange[0]
        xmax=xrange[1]
        yrange=LineageFunctions.getFieldMinMax(points_layer_lin,"lat")
        ymin=yrange[0]
        ymax=yrange[1]
        env.extent = arcpy.Extent(xmin,ymin,xmax,ymax)
        
        # export a kml of the sites        
        group_nospaces = string.replace(group," ","_")
        export_lyr = points_layer_lin
        export_lyr.name = group_nospaces + "_lineages"
        export_lyr.showLabels = True
        lblClass = export_lyr.labelClasses[0]
        lblClass.expression="[Lineage]"
        kmz_name = export_lyr.name + ".kmz"
        arcpy.LayerToKML_conversion(points_layer_lin, kmz_name, layer_output_scale=1, is_composite='NO_COMPOSITE', boundary_box_extent= env.extent)
        
    print "\n   ****************\n   * FINISHED lineage export to kmz for:",group ,"*\n   ****************"
    
# NOW REPEAT THE PROCESS FOR ALL SITES
print "\nAbout to export species sites.\n"

with open(species_site_filename, 'rb') as csvfile:
    species_csv = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_MINIMAL)
    rownum = 0
    GroupList=[]
    Lat=[]
    Long=[]
    
    for row in species_csv:
        # Skip header row.
        if rownum == 0:
            rownum += 1
        else:
            if row[0] not in GroupList and len(row[0])>0:
                GroupList.append(row[0])
    # so now we have a list of model groups, excluding rows with null coordinates

# Make XY Event Layer
try:
    points_layer_sp = "sites"
    layer_result = arcpy.MakeXYEventLayer_management(species_site_filename, "long", "lat", points_layer_sp, spRef)
    points_layer_sp = layer_result[0]
except:
   # If an error occurred print the message to the screen
    print arcpy.GetMessages()
    
# export the layer as a feature class
outLocation = env.workspace
outFeatureClass_sp = genus + "_all_points"

try:
    arcpy.env.overwriteOutput=True
    arcpy.CopyFeatures_management(points_layer_sp, outFeatureClass_sp, "", "0", "0", "0")
except:
    print "\n" + outFeatureClass + " could not be created, or already exists."
points_fc = outLocation + "/" + outFeatureClass_sp

# Loop through the list of Model Groups
for group in GroupList:
    if group != 0 and group != "":
    #if group == "Saproscincus challengeri":   # this is a temp line to model just one group!!
    
        print "\nStarting group " + group + "\n"
        
        groupDefQuery = "[model_group] = '" + group + "'"
        points_layer_sp.definitionQuery = groupDefQuery
               
        # get the extent of the points for the model group
        xrange=LineageFunctions.getFieldMinMax(points_layer_sp,"long")
        xmin=xrange[0]
        xmax=xrange[1]
        yrange=LineageFunctions.getFieldMinMax(points_layer_sp,"lat")
        ymin=yrange[0]
        ymax=yrange[1]
        env.extent = arcpy.Extent(xmin,ymin,xmax,ymax)
        
        # export a kml of the sites        
        group_nospaces = string.replace(group," ","_")
        export_lyr = points_layer_sp
        export_lyr.name = group_nospaces
        export_lyr.showLabels = False
        #lblClass = export_lyr.labelClasses[0]
        #lblClass.expression="[model_group]"
        kmz_name = group_nospaces + "_sites.kmz"
        arcpy.LayerToKML_conversion(points_layer_sp, kmz_name, layer_output_scale=1, is_composite='NO_COMPOSITE', boundary_box_extent= env.extent)
        
    print "\n   ****************\n   * FINISHED species site export to kmz for:",group ,"*\n   ****************"

print "\nDeleting temporary shapefiles\n"
arcpy.Delete_management(outFeatureClass + ".shp")
arcpy.Delete_management(outFeatureClass_sp + ".shp")

print"Finished."