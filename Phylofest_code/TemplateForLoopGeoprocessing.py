import arcpy, sys, os
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

### PARAMETERS ###

points = "C:\\Users\\u3579238\\GISData\\Helping\\Sally\\WyuldaPoints.shp"  # set the point shapefile
buffer_dist = 3   # the buffer distance
source_location = "c:\\Users\\u3579238\\GISData\\EnvironmentGrids\\AusGDMGrids\\maskedgrids\\" # where the original environment grids are stored
target_location = "c:\\Users\\u3579238\\GISData\\Helping\\Sally\\Env_grids_clipped\\" # where the clipped grids go

# create the buffer grid and set as a mask
env.workspace = target_location
env.snapRaster= source_location + "bio01"

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

    #check if the name is a layer to use - in this case bio01 - bio19
    if not ((str.isdigit(str(dataset[-2:])) and int(dataset[-2:]) > 19) or dataset == "bio15_580"):
                    
        # finally, extract the required part of the grid, and write it to specified folder, with a suffix
        maskedgrid = arcpy.sa.ExtractByMask(dataset,pointbuffer)
        new_grid_name = dataset + "_msk.asc"
        print dataset, new_grid_name
        arcpy.RasterToASCII_conversion(maskedgrid, target_location + new_grid_name)

