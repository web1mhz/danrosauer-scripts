import arcpy, sys, numpy, csv, os, string, math
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

# extent limits (to deal with different extents between grids)
Australia_extent = (112.9,-43.75,153.64,-9)

buffer_dist = 2.5   # the buffer distance
source_location = "c:\\Users\\u3579238\\GISData\\EnvironmentGrids\\AusGDMGrids\\maskedgrids\\" # where the original environment grids are stored
target_location = "species_models\\clipped_grids\\" + genus + "\\" # where the clipped grids go

# bioclim layers 1 to 19 are matched automatically.  Any other layers to use are listed here.
extra_layers = ["twi3se_01deg","clay30e_01deg","slope","geollmeanage"]
##################

print "Starting to create clipped enviornment grids with a " + str(buffer_dist) + " buffer."
print "Target: " + target_location

# create the output folder
target_location = base_dir + target_location
if not os.path.exists(target_location):
    os.makedirs(target_location)

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
                        pass
                except:
                    pass
    
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
                    pass
            except:
                pass

#and combine them
sites = species_sites + sequence_sites

os.chdir(target_location)
csv.register_dialect("myDialect",lineterminator="\n")
combined_sites_csv = base_dir + combined_sites_folder + combined_sites_csv
with open(combined_sites_csv,"w") as f:
    new_csv = csv.writer(f, delimiter=',',dialect="myDialect")
    myheader = ("model_group","lat","long")
    new_csv.writerow(myheader)
    new_csv.writerows(sites)
del species_sites, sequence_sites, f

# Make XY Event Layer
try:
    sites_path = string.replace(combined_sites_csv,"\\","/")
    spRef = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
    points_layer = "species_sites"
    arcpy.MakeXYEventLayer_management(sites_path, "long", "lat", points_layer, spRef)
    #arcpy.SaveToLayerFile_management(points_layer, "sites", "ABSOLUTE")    
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
#arcpy.env.extent=arcpy.Extent(pointbuffer)
datasets = arcpy.ListRasters("*","GRID")

for dataset in datasets:

    #check if the name is a layer to use - in this case bio01 - bio19
    if (dataset in extra_layers) or (dataset[0:3] == "bio" and (str.isdigit(str(dataset[-2:])) and int(dataset[-2:]) <= 19)):
                    
        # finally, extract the required part of the grid, and write it to specified folder, with a suffix
        maskedgrid = arcpy.sa.ExtractByMask(dataset,pointbuffer)
        new_grid_name = dataset + "_msk.asc"
        print dataset, new_grid_name
        arcpy.RasterToASCII_conversion(maskedgrid, target_location + new_grid_name)
        
# delete temporary files
grid_path = string.replace(target_location,"\\","/") + "point_buf"
arcpy.Delete_management(grid_path)

print "\nFinished creating buffered grids in " + target_location
