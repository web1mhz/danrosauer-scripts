import arcpy, sys, numpy, csv, os
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

### PARAMETERS ###
os.linesep ="\n"

base_dir = "c:\\Users\\u3579238\\Work\Phylofest\\Models\\skinks\\"
species_site_filename = "species_sites\\Saproscincus_ALA.csv"
sequence_site_filename = "sequence_sites\\Saproscincus_lin_loc.csv"
points = "Diporiphora_ALA.shp"  # set the point shapefile
buffer_dist = 4   # the buffer distance
source_location = "c:\\Users\\u3579238\\GISData\\EnvironmentGrids\\AusGDMGrids\\maskedgrids\\" # where the original environment grids are stored
target_location = "species_models\\clipped_grids\\Saproscincus_gc\\" # where the clipped grids go
##################

# create the buffer grid and set as a mask
target_location = base_dir + target_location
if not os.path.exists(target_location):
    os.makedirs(target_location)

env.workspace = target_location
env.snapRaster= source_location + "bio1"

#load the species site coordinates
species_site_filename = base_dir + species_site_filename
species_sites = numpy.genfromtxt(species_site_filename, delimiter=',',usecols = (1,2),names=True)
#load the sequenced site coordinates
sequence_site_filename = base_dir + sequence_site_filename
sequenced_sites = numpy.genfromtxt(sequence_site_filename, delimiter=',',usecols = (5,6),names=True)
#and combine them
sites = numpy.append(species_sites,sequenced_sites,0)
sites = sites[(sites[:,1] > 100) & (sites[:,0] > -45)]  #removes null coordinates and some wrong ones

os.chdir(target_location)
csv.register_dialect("myDialect",lineterminator="\n")
new_csv = csv.writer(open("sites_temp.csv","w"), delimiter=',',dialect="myDialect")
myheader = ("lat","long")
new_csv.writerow(myheader)
new_csv.writerows(list(sites))
del species_sites
del sequenced_sites

# Make XY Event Layer
sites_path = target_location + "sites_temp.csv"
spRef = "C:\Program Files (x86)\ArcGIS\Desktop10.0\Coordinate Systems\Geographic Coordinate Systems\World\WGS 1984.prj"
arcpy.MakeXYEventLayer_management(sites_path, "long", "lat", "sites.shp", spRef)

# Process: Copy Features
arcpy.CopyFeatures_management(table_Layer2, point3_shp, "", "0", "0", "0")

#get the points layer extent, and calculate the buffer extent
points_properties = arcpy.Describe(points)
points_extent = points_properties.extent
extent_buffer = buffer_dist * 1.25
new_extent = [points_extent.xmin - extent_buffer, points_extent.ymin - extent_buffer, points_extent.xmax + extent_buffer, points_extent.ymax + extent_buffer]
env.extent= arcpy.Extent(new_extent[0],new_extent[1],new_extent[2],new_extent[3]) # set the analysis extent

pointbuffer = arcpy.sa.EucDistance(points,buffer_dist,0.01)
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

