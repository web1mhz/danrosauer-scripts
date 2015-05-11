import arcpy
import sys, numpy, csv, os, string, math, subprocess, fnmatch
from datetime import datetime
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

### PARAMETERS ###
os.linesep ="\n"

#higher_taxon = "skinks"
genus_list = ["Gehyra"]  # genus could refer to any group being handled as a set

# a changeable list to allow for species in the dataset to be skipped
named_species   = ["Gehyra_multiporosa_grp","Gehyra_nana","Gehyra_sp_nov"]
use_list        = "do"    #specify whether to:
                            #do only the named species (use_list="do")
                            #skip the named species (use_list="skip")
                            #do all the species in the data and ignore the list (use_list="" or anything else);

#base_dir = "C:\\Users\\u3579238\\work\\Phylofest\\Models\\" + higher_taxon + "\\"
base_dir = "C:\\Users\\u3579238\\work\\AMT\\Models\\"
combined_sites_folder    = "species_sites\\"
maxent_loc = "C:\\Users\\u3579238\\Work\\Phylofest\Models\\maxent.jar"

Australia_extent = (112.9,-43.75,153.64,-9)

buffer_dist = 3   # the buffer distance
source_location = "c:\\Users\\u3579238\\GISData\\EnvironmentGrids\\AusGDMGrids\\maskedgrids\\" # where the original environment grids are stored
extra_layers = ["twi3se_01deg","slope","geollmeanage"]

use_bias_grid  = False
bias_grid_name = "geckoes_bias_grid"
#bias_grid_name = "skink_samplesites" 
bias_grid_source_loc = base_dir + "species_models\\bias_files\\bias.gdb\\"
bias_grid_temp_loc   = base_dir + "species_models\\bias_files\\"
maxent_replicates = 25
jackknife = False
processor_threads = 8  # this sets how many processors MaxEnt can use.

output_gdb_name = "maxent_models.gdb" 
model_suffix = "_median"

done_list = []

for genus in genus_list:
    
    #species_site_filename = "species_sites\\" + genus + "_ALA.csv"
    species_site_filename = "\\species_sites\\" + genus + "_species_sites.csv"
    combined_sites_csv = genus + "_maxent.csv" #name of the new csv file of species,lat,long to be created
    lineage_site_filename = "lineage_sites\\" + genus + "_lin_loc.csv"
    
    maxent_model_base = base_dir + "species_models\\maxent\\" + genus
    
    target_location = "species_models\\clipped_grids\\" + genus + "\\" # where the clipped grids go
    
    # bioclim layers 1 to 19 are matched automatically.  Any other layers to use are listed here.
        
    ##################
        
    print "Starting to create clipped environment grids with a buffer of " + str(buffer_dist)
    print "Target: " + target_location
    print datetime.now()
    
    # create the output folder for clipped grids
    target_location = base_dir + target_location
    if not os.path.exists(target_location):
        os.makedirs(target_location)
        
    # create the output geodatabase for species models, if needed
    if not os.path.exists(maxent_model_base):
        os.makedirs(maxent_model_base)
    arcpy.env.overwriteOutput=False
    try:
        maxent_model_GDB_path = base_dir + "species_models\\maxent\\"
        maxent_model_GDB_path = string.replace(maxent_model_GDB_path,"\\","/")
        arcpy.CreateFileGDB_management(maxent_model_GDB_path, output_gdb_name)
    except:
        print "\nFile geodatabase "+maxent_model_GDB_path + output_gdb_name + " exists, or could not be created.\n"
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
        usecol = -1
    
        for row in sequence_csv:
            # Save header row.
            if rownum == 0:
                header = row
                rownum += 1
                columns = range(len(header))            
    
                # find the column called 'use'
                for k in columns:
                    if str.lower(header[k]) == "use":
                        usecol = k
                        break
                # find the column called 'Latitude'
                for k in columns:
                    if str.lower(header[k]) == "latitude" or str.lower(header[k]) == "lat":
                        lat_col = k
                        break
                # find the column called 'Longitude'
                for k in columns:
                    if str.lower(header[k]) == "longitude"or str.lower(header[k]) == "long":
                        long_col = k
                        break
                # find the column called 'Species'
                for k in columns:
                    if str.lower(header[k]) == "species":
                        species_col = k
                        break                  
            else:
                rownum += 1
                if (row[usecol] == '1'): # only proceed where 'use' = 1 and x and y are within defined limits
                    try:        # if the lat or long can't be converted to a number, then skip that row, by not incremeting rownum
                        num = (float(row[lat_col]))
                        try:
                            # code gets to here for valid lat and long, so other steps can go here too
                            num = (float(row[long_col]))
                            row[species_col] = string.replace(row[species_col]," ","_")  #replace spaces in taxon name with _
                            row[lat_col] = float(row[lat_col])  #change x,y coords from text to numbeLatituder
                            row[long_col] = float(row[long_col])
                            if (row[lat_col] > ylim_min_points) and (row[long_col] > xlim_min_points): # add row if x and y are within defined limits
                                species_sites.append(row[:3])
                        except:
                            pass
                    except:
                        pass
        
    #load the sequenced site coordinates
    lineage_site_filename = base_dir + lineage_site_filename
    with open(lineage_site_filename, 'rb') as csvfile:
        sequence_csv = csv.reader(csvfile, delimiter=',')
        rownum = 0
        usecol = 0
        sequence_sites = []
    
        for row in sequence_csv:
            # Save header row.
            if rownum == 0:
                header = row
                rownum += 1
                
                # find the row called 'use'
                columns = range(len(header))
                for k in columns:
                    if str.lower(header[k]) == "use":
                        usecol = k
                        break
            else:
                if (usecol == 0 or row[usecol] == '1'): # if there is a use column, only proceed where 'use' = 1 and x and y are within defined limits
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
                            pass
                    except:
                        pass
    
    #and combine them
    #sites = species_sites + sequence_sites
    sites = species_sites
    
    #list the model_groups and remove spaces from taxon names
    rownum = 0
    model_groups = []
    for row in sites:
        row[0] = string.replace(row[0]," ","_") #replace space in taxon name with _
        rownum += 1
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
            new_extent = [points_extent.XMin - extent_buffer, points_extent.YMin - extent_buffer, points_extent.XMax + extent_buffer, points_extent.YMax + extent_buffer]
        
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
    
            # include the bias grid in rasters to be clipped, copied
            if use_bias_grid:
                
                # if the bias grid is ascii, convert it to raster
                splitname = string.rsplit(bias_grid_name,".",1)
                if len(splitname) > 1:
                    if splitname[1] == "asc":
                        new_bias_grid_name = "bias_" + model_group[:4]
                        arcpy.ASCIIToRaster_conversion(bias_grid_temp_loc + bias_grid_name, data_type="INTEGER", out_raster=bias_grid_temp_loc + new_bias_grid_name)
                        bias_grid_name = new_bias_grid_name
                        del(new_bias_grid_name)
                
                env.workspace = bias_grid_source_loc
                env.cellSize=(bias_grid_name)
                
                ## combine the points used in the model, to ensure all are in the grid
                #points_layer_grid = arcpy.PointToRaster_conversion(in_features=points_layer, value_field="model_group")
                ##new_bias_grid = arcpy.gp.RasterCalculator_sa("Con(Con(IsNull("extra_samplesites_test"),0,1) + Con(IsNull("skink_samplesites"),0,1)>0,1)""","C:/Users/u3579238/Work/Phylofest/Models/skinks/species_models/bias_files/bias.gdb/skink_samplesites_fulltest")
                #new_bias_grid = arcpy.sa.Con(arcpy.sa.Con(arcpy.sa.IsNull(points_layer_grid),0,1) + arcpy.sa.Con(arcpy.sa.IsNull(bias_grid_name),0,1)>0,1)
                new_bias_grid = bias_grid_name
                
                #extract the required part of the grid, and write it to specified folder, with a suffix
                maskedgrid = arcpy.sa.ExtractByMask((new_bias_grid),pointbuffer)
                masked_bias_grid_name = bias_grid_name + "_" + model_group_sp + "_msk.asc"
                print "\nBias grid for " + model_group + ": ", new_grid_name
                arcpy.RasterToASCII_conversion(maskedgrid, bias_grid_temp_loc + masked_bias_grid_name)
                    
            #RUN THE MAXENT MODEL HERE
            # create the output folder
            if not os.path.exists(maxent_model_base):
                os.makedirs(maxent_model_base)
                
            maxent_call = "java -mx1024m -cp " + maxent_loc + " density.MaxEnt nowarnings noprefixes novisible outputdirectory=" + maxent_model_base +  " samplesfile=" + model_group_sites_csv + " environmentallayers=" + target_location + " replicates=" + str(maxent_replicates) + " autorun randomseed threads=" + str(processor_threads)
            
            if jackknife:
                maxent_call += " jackknife "
            else:
                maxent_call += " nojackknife "
            
            if use_bias_grid:
                maxent_call += " biasfile=" + bias_grid_temp_loc + masked_bias_grid_name + " biastype=3"
            
            print maxent_call  # drop this line after debug complete
            
            print "\nAbout to start maxent model for: " + model_group + "  replicates: " + str(maxent_replicates) + "\nat" + str(datetime.now()) + "\n"
            
            try:
                subprocess.call(maxent_call)
             
                #COPY THE RESULT TO GDB
                # set the geoprocessing environment
                model_gdb = maxent_model_GDB_path + output_gdb_name
                env.workspace  = model_gdb
            
                # import the maxent model result from ascii to ESRI gdb
                maxent_model = maxent_model_base + "\\" + string.replace(model_group," ","_") + model_suffix + ".asc"
                out_raster = model_gdb+"/" + string.replace(model_group," ","_")
                arcpy.ASCIIToRaster_conversion(maxent_model, out_raster, "FLOAT")
                print "\nImported " + string.replace(model_group," ","_") + ".asc" + " to " + out_raster
            
                print "\nFinished model import to GDB for: " + model_group + "\n"
            except:
                print "\n*******************\nModel failed for ", model_group, "\n\n"
                
            # delete ArcGIS temporary files
            print "Now deleting buffered grids for: " + model_group
            env.workspace = target_location
            wildcard = genus + "*"
            keep_datasets = arcpy.ListRasters(wildcard)
            all_datasets =  arcpy.ListRasters("*")
            for dataset in all_datasets:
                if dataset not in keep_datasets:
                    arcpy.Delete_management(dataset)
                    
            # delete temporary buffered bias grid
            if use_bias_grid:
                print "Deleting clipped bias grid: " + masked_bias_grid_name
                arcpy.Delete_management(bias_grid_temp_loc + masked_bias_grid_name)
                    
            # delete Maxent temporary files
            print "Now deleting redundant Maxent files for: " + model_group
            os.chdir(maxent_model_base)
            files = os.listdir(maxent_model_base)
            files_to_delete = fnmatch.filter(files,'*_[0-9]*')
            del_count = 0
            for file in files_to_delete:
                os.remove(file)
                del_count +=1
    
            plots_dir = maxent_model_base + "\\plots"
            os.chdir(plots_dir)
            files = os.listdir(plots_dir)
            files_to_delete = fnmatch.filter(files,'*_[0-9]*')
            for file in files_to_delete:
                os.remove(file)
                del_count +=1
            
            print "Deleted " + str(del_count) + " files created by Maxent\n"
            print "\n*********************************\nFinished models for " + model_group +"\nat " + str(datetime.now()) + "\n*********************************\n"
            done_item = (genus + "  " + model_group)
            done_list.append(done_item)
        else:
            print "\n****************************\nSkipping " + model_group + "\n****************************\n"
            
    print "\n****************************\nFinished all models for " + genus + " at " + str(datetime.now()) + "\n****************************\n"

print "\n****************************\nFinished ALL models at " + str(datetime.now()) + "\n****************************\n"
for item in done_list:
    print item
print "\n****************************\nFinished ALL models at " + str(datetime.now()) + "\n****************************\n"
