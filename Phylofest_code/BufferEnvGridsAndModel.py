import arcpy, sys, numpy, csv, os, string, math, subprocess
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

### PARAMETERS ###
os.linesep ="\n"

genus = "Lampropholis"  # genus could refer to any group being handled as a set
higher_taxon = "skinks"
base_dir = "C:\\Users\\u3579238\\work\\Phylofest\\Models\\" + higher_taxon + "\\"
species_site_filename = "species_sites\\" + genus + "_ALA.csv"
combined_sites_folder    = "species_sites\\"
combined_sites_csv = genus + "_maxent.csv" #name of the new csv file of species,lat,long to be created
sequence_site_filename = "sequence_sites\\" + genus + "_lin_loc.csv"
maxent_loc = "C:\\Users\\u3579238\\Work\\Phylofest\Models\\maxent.jar"
maxent_model_base = base_dir + "species_models\\maxent\\" + genus + "\\"

# extent limits (to deal with different extents between grids)
Australia_extent = (112.9,-43.75,153.64,-9)

buffer_dist = 2.5   # the buffer distance
source_location = "c:\\Users\\u3579238\\GISData\\EnvironmentGrids\\AusGDMGrids\\maskedgrids\\" # where the original environment grids are stored
target_location = "species_models\\clipped_grids\\" + genus + "\\" # where the clipped grids go

# bioclim layers 1 to 19 are matched automatically.  Any other layers to use are listed here.
extra_layers = ["twi3se_01deg","clay30e_01deg","slope","geollmeanage"]

output_gdb_name = "maxent_models.gdb"
model_suffix = "_median"
maxent_replicates = 25

# a changeable list to allow for species in the dataset to be skipped
named_species   = ["Lampropholis_amicula", "Lampropholis_robertsi", "Lampropholis_guichenoti"]
use_list        = "do"  #specify whether to:
                            #do the named species (use_list="do")
                            #skip the named species (use_list="skip")
                            #do all the species in the data and ignore the list (use_list="" or anything else);

##################

print "Starting to create clipped enviornment grids with a " + str(buffer_dist) + " buffer."
print "Target: " + target_location

# create the output folder for clipped grids
target_location = base_dir + target_location
if not os.path.exists(target_location):
    os.makedirs(target_location)
    
# create the output geodatabase for species models, if needed
if not os.path.exists(maxent_model_base):
    os.makedirs(maxent_model_base)
arcpy.env.overwriteOutput=False
try:
    maxent_model_base_ESRI = string.replace(maxent_model_base,"\\","/")
    arcpy.CreateFileGDB_management(maxent_model_base_ESRI, output_gdb_name)
except:
    print "\nFile geodatabase "+maxent_model_base_ESRI+ " " + output_gdb_name + " exists, or could not be created.\n"
arcpy.env.overwriteOutput=True    

env.workspace = target_location
env.snapRaster= source_location + "bio01"

# coordinate limits for filtering points
xlim_min_points = Australia_extent[0]
ylim_min_points = Australia_extent[1]

#load the species site coordinates
species_site_filename = base_dir + species_site_filename
with open(species_site_filename, 'rb') as csvfile:
    sequence_csv = csv.reader(csvfile, delimiter=',')
    rownum = 0
    species_sites = []
    num=0

    for row in sequence_csv:
        # Save header row.
        if rownum == 0:
            header = row
            rownum += 1
        else:
            if (row[4] == '1'): # only proceed where 'use' = 1 and x and y are within defined limits
                try:        # if the lat or long can't be converted to a number, then skip that row, by not incremeting rownum
                    num = (float(row[1]))
                    try:
                        # code gets to here for valid lat and long, so other steps can go here too
                        num = (float(row[2]))
                        row[0] = string.replace(row[0]," ","_")  #replace spaces in taxon name with _
                        row[1] = float(row[1])  #change x,y coords from text to number
                        row[2] = float(row[2])
                        if (row[1] > ylim_min_points) and (row[2] > xlim_min_points): # add row if x and y are within defined limits
                            species_sites.append(row[:3])
                    except:
                        1==1
                except:
                    1==1
    
#load the sequenced site coordinates
sequence_site_filename = base_dir + sequence_site_filename
with open(sequence_site_filename, 'rb') as csvfile:
    sequence_csv = csv.reader(csvfile, delimiter=',')
    rownum = 0
    sequence_sites = []

    for row in sequence_csv:
        # Save header row.
        if rownum == 0:
            header = row
            rownum += 1
        else:
            try:        # if the lat or long can't be converted to a number, then skip that row, by not incremeting rownum
                num = (float(row[5]))
                try:
                    # code gets to here for valid lat and long, so other steps can go here too
                    num = (float(row[6]))
                    row[5] = float(row[5])  #change x,y coords from text to number
                    row[6] = float(row[6])
                    if (row[5] > ylim_min_points) and (row[6] > xlim_min_points): # add row if x and y are within defined limits
                        sequence_sites.append([row[3],row[5],row[6]])
                    rownum += 1
                except:
                    1==1
            except:
                1==1

#and combine them
sites = species_sites + sequence_sites

#list the model_groups and remove spaces from taxon names
model_groups = []
for row in sites:
    row[0] = string.replace(row[0]," ","_") #replace space in taxon name with _
    if row[0] not in model_groups:
        model_groups.append(row[0])
        
    #if needed, create a folder for location files in the genus
    sites_dir = base_dir + combined_sites_folder + genus
    if not os.path.exists(sites_dir):
        os.makedirs(sites_dir)
        
# write the combined location record list to file
os.chdir(target_location)
csv.register_dialect("myDialect",lineterminator="\n")
combined_sites_csv = base_dir + combined_sites_folder + combined_sites_csv
with open(combined_sites_csv,"w") as f:
    new_csv = csv.writer(f, delimiter=',',dialect="myDialect")
    myheader = ("model_group","lat","long")
    new_csv.writerow(myheader)
    new_csv.writerows(sites)
    
del species_sites, sequence_sites, f

# restrict the model_groups to particular species based on the names_species parameter
if use_list == "do":
    model_groups = list(set(named_species).intersection(set(model_groups)))
elif use_list == "skip":
    model_groups = list(set(model_groups).difference(set(named_species)))

# loop through the model_groups
for model_group in model_groups:
    if model_group.find("_") >= 0:   # only model where the taxon name is binomial
        group_sites=[]
        for row in sites:
            if row[0] == model_group:
                group_sites.append(row)
                
        # write the model_group location record list to file
        model_group_sites_csv = base_dir + combined_sites_folder + genus + "\\" + model_group + ".csv"
        with open(model_group_sites_csv,"w") as f:
            new_csv = csv.writer(f, delimiter=',',dialect="myDialect")
            myheader = ("model_group","lat","long")
            new_csv.writerow(myheader)
            new_csv.writerows(group_sites)
    
        # Make XY Event Layer for model group
        try:
            sites_path = string.replace(model_group_sites_csv,"\\","/")
            spRef = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
            points_layer = "species_sites"
            arcpy.MakeXYEventLayer_management(sites_path, "long", "lat", points_layer, spRef)
        except:
           # If an error occurred print the message to the screen
            print arcpy.GetMessages() 
    
        #get the points layer extent, and calculate the buffer extent
        points_properties = arcpy.Describe(points_layer)
        points_extent = points_properties.extent
        extent_buffer = buffer_dist * 1.25
        new_extent = [points_extent.xmin - extent_buffer, points_extent.ymin - extent_buffer, points_extent.xmax + extent_buffer, points_extent.ymax + extent_buffer]
    
        # but where the extended buffer goes beyond the extent of Australia, limit to Australia
        xmin=max(round(new_extent[0],2),Australia_extent[0])
        ymin=max(round(new_extent[1],2),Australia_extent[1])
        xmax=min(round(new_extent[2],2),Australia_extent[2])
        ymax=min(round(new_extent[3],2),Australia_extent[3])
        env.extent = arcpy.Extent(xmin,ymin,xmax,ymax)
    
        pointbuffer = arcpy.sa.EucDistance(points_layer,buffer_dist,0.01)
        pointbuffer.save("point_buf")
    
        env.workspace = source_location
        datasets = arcpy.ListRasters("*","GRID")
    
        for dataset in datasets:
        
            #check if the name is a layer to use - in this case bio01 - bio19
            if (dataset in extra_layers) or (dataset[0:3] == "bio" and (str.isdigit(str(dataset[-2:])) and int(dataset[-2:]) <= 19)):
                            
                # finally, extract the required part of the grid, and write it to specified folder, with a suffix
                maskedgrid = arcpy.sa.ExtractByMask(dataset,pointbuffer)
                model_group_sp = model_group[(model_group.find("_")+1):] # use just the species name in the grid name (not the whole binomial)
                new_grid_name = dataset + "_" + model_group_sp + "_msk.asc"
                print model_group, dataset, new_grid_name
                arcpy.RasterToASCII_conversion(maskedgrid, target_location + new_grid_name)
                
        #RUN THE MAXENT MODEL HERE
        # create the output folder
        if not os.path.exists(maxent_model_base):
            os.makedirs(maxent_model_base)
        maxent_call = "java -mx1024m -cp " + maxent_loc + " density.MaxEnt nowarnings noprefixes novisible jackknife outputdirectory=" + maxent_model_base +  " samplesfile=" + model_group_sites_csv + " environmentallayers=" + target_location + " replicates=" + str(maxent_replicates) + " autorun randomseed"
        
        print "\nAbout to start maxent model for: " + model_group + "  replicates: " + str(maxent_replicates) + "\n"
        subprocess.call(maxent_call)
          
        #COPY THE RESULT TO GDB
        # set the geoprocessing environment
        model_gdb = maxent_model_base_ESRI + output_gdb_name
        env.workspace  = model_gdb
    
        # import the maxent model result from ascii to ESRI gdb
        maxent_model = maxent_model_base + "\\" + string.replace(model_group," ","_") + model_suffix + ".asc"
        out_raster = model_gdb+"/" + string.replace(model_group," ","_")
        arcpy.ASCIIToRaster_conversion(maxent_model, out_raster, "FLOAT")
        print "\nImported " + string.replace(model_group," ","_") + ".asc" + " to " + out_raster
    
        print "\nFinished model import to GDB for: " + model_group + "\n"
        print "Now deleting buffered grids for: " + model_group
        
        # delete temporary files
        env.workspace = target_location
        wildcard = genus + "*"
        keep_datasets = arcpy.ListRasters(wildcard)
        all_datasets =  arcpy.ListRasters("*")
        for dataset in all_datasets:
            if dataset not in keep_datasets:
                arcpy.Delete_management(dataset)
        
        print "\n****************************\nFinished models for " + model_group + "\n****************************\n"
        
    else:
        print "\n****************************\nSkipping " + model_group + "\n****************************\n"
        
print "\n****************************\nFinished all models for " + genus + "\n****************************\n"

