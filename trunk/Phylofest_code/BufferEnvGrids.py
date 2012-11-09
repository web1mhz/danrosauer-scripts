import arcpy, sys, numpy, csv, os, string
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

### PARAMETERS ###
os.linesep ="\n"

base_dir = "C:\\Users\\u3579238\\work\\Phylofest\\Models\\skinks\\"
species_site_filename = "species_sites\\Saproscincus_ALA.csv"
sequence_site_filename = "sequence_sites\\Saproscincus_lin_loc.csv"
buffer_dist = 3.5   # the buffer distance
source_location = "c:\\Users\\u3579238\\GISData\\EnvironmentGrids\\AusGDMGrids\\maskedgrids\\" # where the original environment grids are stored
target_location = "species_models\\clipped_grids\\Saproscincus\\" # where the clipped grids go
##################

# create the output folder
target_location = base_dir + target_location
if not os.path.exists(target_location):
    os.makedirs(target_location)

env.workspace = target_location
env.snapRaster= source_location + "bio1"

#load the species site coordinates
species_site_filename = base_dir + species_site_filename

try:     #removes coordinates where the 5th column (col num = 4) is not set to 1
    species_sites = numpy.genfromtxt(species_site_filename, delimiter=',',usecols = (1,2,4),names=True)
    species_sites = species_sites[(species_sites[:,2] == 1),:2]  #removes coordinates where use is not set to 1
except:  # or if there is no 5th column, just loads 1,2 for x,y
    species_sites = numpy.genfromtxt(species_site_filename, delimiter=',',usecols = (1,2),names=True)
    
#load the sequenced site coordinates
sequence_site_filename = base_dir + sequence_site_filename
sequenced_sites = numpy.genfromtxt(sequence_site_filename, delimiter=',',usecols = (5,6),names=True)
#and combine them
sites = numpy.append(species_sites,sequenced_sites,0)
sites = sites[(sites[:,1] > 100) & (sites[:,0] > -45)]  #removes null coordinates and some wrong ones

os.chdir(target_location)
csv.register_dialect("myDialect",lineterminator="\n")
sites_temp_csv = "sites_temp.csv"
with open(sites_temp_csv,"w") as f:
    new_csv = csv.writer(f, delimiter=',',dialect="myDialect")
    myheader = ("lat","long")
    new_csv.writerow(myheader)
    new_csv.writerows(list(sites))
del species_sites
del sequenced_sites
del f

# Make XY Event Layer
try:
    sites_path = string.replace(target_location,"\\","/") + sites_temp_csv
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
env.extent= arcpy.Extent(new_extent[0],new_extent[1],new_extent[2],new_extent[3]) # set the analysis extent

pointbuffer = arcpy.sa.EucDistance(points_layer,buffer_dist,0.01)
pointbuffer.save("point_buf")

env.workspace = source_location
arcpy.env.extent=arcpy.Extent(pointbuffer)
datasets = arcpy.ListRasters("*bio*","GRID")

for dataset in datasets:

    #check if the name is a layer to use - in this case bio1 - bio19
    if not ((str.isdigit(str(dataset[-2:])) and int(dataset[-2:]) > 19) or dataset == "bio15_580"):
                    
        # finally, extract the required part of the grid, and write it to specified folder, with a suffix
        maskedgrid = arcpy.sa.ExtractByMask(dataset,pointbuffer)
        new_grid_name = dataset + "_msk.asc"
        print dataset, new_grid_name
        arcpy.RasterToASCII_conversion(maskedgrid, target_location + new_grid_name)
        
# delete temporary files
#os.remove(sites_path)
grid_path = string.replace(target_location,"\\","/") + "point_buf"
arcpy.Delete_management(grid_path)

print "\nFinished creating buffered grids in " + target_location
