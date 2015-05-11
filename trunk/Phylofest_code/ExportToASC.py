# Dan Rosauer November 2012

# reads in the lineage distributions from an ESRI file GDB and exports an .asc for each lineage.

import arcpy, sys, os, string
from arcpy import env
import arcpy.sa
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

sys.path.append("C:\\Users\\Dan\Work\AMT\\")
import LineageFunctions

### PARAMETERS ###
higher_taxon = "geckoes"
#genus_list = ["Carlia","Eulamprus","Glaphyromorphus","Gnypetoscincus","Saproscincus","Lampropholis"]
genus_list = ["Heteronotia"]

base_dir = "C:\\Users\\u3579238\\work\\AMT\\Models\\"
target_location = "C:\\Users\\u3579238\\Work\\AMT\\Models\\lineage_models\\asc_new_method\\"
    
for genus in genus_list:
    print "\n" + genus + "\n"

    #lineage_model_dir = base_dir + "lineage_models\\" + genus + "\\results.gdb"
    lineage_model_dir = base_dir + "lineage_models\\lineage_models.gdb"

    arcpy.env.workspace = lineage_model_dir
    arcpy.env.overwriteOutput=True
    gridlist= arcpy.ListRasters()

    coordinateSystem = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"

    for grid in gridlist:
        out_filename = target_location + grid + ".asc"
        arcpy.RasterToASCII_conversion(grid, out_filename)
        arcpy.DefineProjection_management(out_filename, coordinateSystem)

        print grid 

print"Finished."
