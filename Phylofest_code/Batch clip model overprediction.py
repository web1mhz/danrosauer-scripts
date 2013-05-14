# Dan Rosauer May 2013

##STEPS TO ENCODE
#
#in the model correction polygon feature class, create polygons to be used in 3 ways (so far)
#
#1) exclude. For the named species, set areas in the polygon to 0
#2) include. For the named species, set areas outside the polygon to 0
#3) fade.    For the named species, areas outside the polygon, up to the specified distance (D) are downweighted as follows.  Beyond this distance, value is 0:
#	
#	[a] euclidean distance outside the polygon (no maximum)
#	[b] use Con to set eucdist > D to D. (so distance plateaus at D)
#	[c] new model grid = model * ((D - [b]) / D)
#
#      So, the model value is:
#	(a) unchanged within the polygon
#	(b) declines from its original value to 0 over a distance D from the polygon
#	(c) is 0 beyond distance D from the polygon
#
#Environments:
#Use original model for extent, snap raster, grid size
#
#The algorithm will loop through rows in the model correction polygon feature class, and apply each rule to the specified model
#
#This will make the edits repeatable, explicit.

import arcpy, sys, os, math, numpy, csv, os, string
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

sys.path.append("C:\\Users\\u3579238\Work\Phylofest\\")
import LineageFunctions

### PARAMETERS ###
rules_gdb = "C:\\Users\\u3579238\\Work\\Phylofest\\Models\\EditModels.gdb"
rules_fc  = "Model_edits"
higher_taxon = "geckoes"
model_dir = "C:\\Users\\u3579238\\work\\Phylofest\\Models\\" + higher_taxon + "\\species_models\\"
output_location = rules_gdb  # where the edited models go.  Could be different to the rules_gdb
scratch_workspace = "C:\\Users\\u3579238\\Work\\Phylofest\\Models\\"  #scratch workspace is used by ArcGIS for temporary files during analysis
scratch_gdb = "C:\\Users\\u3579238\\Work\\Phylofest\\Models\\scratch.gdb"
grid_resolution = 0.01
 
### create the scratch geodatabase if needed (to store temporary layers during analysis)
#try:
#    arcpy.CreateFileGDB_management(scratch_workspace, scratch_gdb)
#except:
#    print "File geodatabase "+ scratch_gdb + " already exists or creation failed"

## create the output geodatabase if needed
#if not os.path.exists(target_location):
#    os.makedirs(target_location)
#try:
#    target_location_ESRI = string.replace(target_location,"\\","/")
#    arcpy.CreateFileGDB_management(target_location_ESRI, "results.gdb")
#except:
#    print "File geodatabase "+target_location + output_gdb_name + " already exists or creation failed"

# Set the geoprocessing environment
env.workspace = scratch_gdb
arcpy.env.grid_resolution = grid_resolution
# Load the edit rule feauture class
print "\nLoading the edit rules\n"
rules_fc_path = rules_gdb + "\\" + rules_fc
rule_cursor   = arcpy.da.SearchCursor(rules_fc_path, ["OBJECTID","SHAPE@","Action","Species","Distance","higher_taxon"])
for row in rule_cursor:
    rule_poly = row[1]
    rule_action = row[2]
    rule_taxon= row[3]
    rule_distance=row[4]
    rule_higher_taxon = row[5]
    
    if rule_higher_taxon == higher_taxon:
    
        model_in = model_dir + "maxent\\" + rule_taxon.split("_")[0] + "\\maxent_models.gdb\\" + rule_taxon
        env.snapRaster = model_in
        env.extent     = model_in
        env.mask       = model_in
        model_out_name = output_location + "\\" + rule_taxon + "_edited"
        
        # check for fade where distance = 0 (would cause division by 0 error)
        if rule_action == "fade" and rule_distance <= 0:
            rule_action = "include"
        
        # implement an action based on the value of Action
        if rule_action == "fade" and rule_distance > 0:
            print "\nAbout to do " + rule_taxon + ". Fade over " + str(round(rule_distance,4)) + " degrees\n"
            eucdist_poly = arcpy.sa.EucDistance(rule_poly, 10, grid_resolution)
            eucdist_poly_con = arcpy.sa.Con(eucdist_poly >= rule_distance, rule_distance, eucdist_poly)
            model_out = model_in * ((rule_distance - eucdist_poly_con) / rule_distance)
            
        elif rule_action == "include":
            print "\nAbout to do " + rule_taxon + ". Include only areas within polygon\n"
            extract_by_poly = arcpy.sa.ExtractByMask(model_in,rule_poly)
            model_out = arcpy.sa.Con(arcpy.sa.IsNull(extract_by_poly), 0, extract_by_poly)
            
        elif rule_action == "exclude":
            print "\nAbout to do " + rule_taxon + ". Exclude areas within polygon\n"
            extract_by_poly = arcpy.sa.ExtractByMask(model_in,rule_poly)
            model_neg = m= arcpy.sa.Con(arcpy.sa.IsNull(extract_by_poly), 0, extract_by_poly)
            model_out = model_in - model_neg
            
        model_out.save(model_out_name)
        model_out.save(model_in + "_edited")
     
# PUT EDITED RASTER BACK IN ORIGINAL LOCATION

# ADD CODE HERE TO CLEAN UP SCRATCH GBD
   
print "\nModel clipping done!\n"
