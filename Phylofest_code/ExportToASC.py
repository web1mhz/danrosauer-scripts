# Dan Rosauer November 2012

# a quickly edited version of the lineage range method, to read in the species and lineage distributions and export a .kml for each species.

import arcpy, sys, os
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

sys.path.append("C:\\Users\\u3579238\Work\Phylofest\\")
import LineageFunctions

### PARAMETERS ###
genus = "Saproscincus_South"  # genus could refer to any group being handled as a set
higher_taxon = "skinks"
base_dir = "C:\\Users\\u3579238\\work\\Phylofest\\Models\\" + higher_taxon + "\\"
lineage_model_dir = base_dir + "lineage_models\\" + genus + "\\results.gdb"
target_location = "C:\\Users\\u3579238\\Work\\Phylofest\\Models\\combined\\lineage_models\\"

arcpy.env.workspace = lineage_model_dir
gridlist= arcpy.ListRasters()

coordinateSystem = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"

for grid in gridlist:
    out_filename = target_location + grid + ".asc"
    arcpy.RasterToASCII_conversion(grid, out_filename)
    arcpy.DefineProjection_management(out_filename, coordinateSystem)
    
    print grid 

print"Finished."