# Dan Rosauer September 2012

import arcpy, sys, os, math, numpy, csv, os, string
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

### PARAMETERS ###
genus = "Lampropholis"  # genus could refer to any group being handled as a set
higher_taxon = "skinks"
base_dir = "C:\\Users\\u3579238\\work\\Phylofest\\Models\\" + higher_taxon + "\\"
maxent_model_base = base_dir + "species_models\\maxent\\" + genus + "\\"
output_gdb_name = "maxent_models.gdb"
sequence_site_filename = base_dir + "sequence_sites\\" + genus + "_lin_loc.csv"
model_suffix = "_median"

##################
# create the output geodatabase if needed
if not os.path.exists(maxent_model_base):
    os.makedirs(maxent_model_base)
try:
    maxent_model_base_ESRI = string.replace(maxent_model_base,"\\","/")
    arcpy.CreateFileGDB_management(maxent_model_base_ESRI, output_gdb_name)
    
except:
    print "File geodatabase "+maxent_model_base_ESRI+ " " + output_gdb_name + " already exists or creation failed"

# Load the sequence site data
print "\nLoading the sequenced sites\n"

# make each column being used, into a list
# and also create lists of unique Analysis Groups and Lineages
with open(sequence_site_filename, 'rb') as csvfile:
    sequence_csv = csv.reader(csvfile, delimiter=',')
    rownum = 0
    Lat=[]
    Long=[]
    GroupList=[]
    
    for row in sequence_csv:
        # Save header row.
        if rownum == 0:
            header = row
            rownum += 1
        else:

            try:        # if the lat or long can't be converted to a number, then skip that row, by not incremeting rownum
                Lat.append(float(row[5]))
                try:
                    # code gets to here for valid lat and long, so other steps can go here too
                    Long.append(float(row[6]))
                    if row[3] not in GroupList:
                        GroupList.append(row[3])
                    rownum += 1
                except:
                    a=0 # just a do nothing line to keep the syntax working
            except:
                a=0 # just a do nothing line to keep the syntax working
                
    # so now we have a list for each column, excluding rows with null coordinates
del Lat
del Long

# set the geoprocessing environment
model_gdb = maxent_model_base_ESRI + output_gdb_name
env.workspace  = model_gdb

# Loop through the list of Model Groups
for group in GroupList:
    
    print "\n" + group

    # import the maxent model result from ascii to ESRI gdb
    maxent_model = maxent_model_base + "\\" + string.replace(group," ","_") + model_suffix + ".asc"
    out_raster = model_gdb+"/" + string.replace(group," ","_")
    arcpy.ASCIIToRaster_conversion(maxent_model, out_raster, "FLOAT")
    print "\nImported " + string.replace(group," ","_") + ".asc" + " to " + out_raster

print "\nFinished model import\n"