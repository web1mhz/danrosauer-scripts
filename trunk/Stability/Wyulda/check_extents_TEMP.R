
rm(list=ls())
library(raster)
library(SDMTools)

################################################################################
#define directories
source.dir = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/mxe/temp'
################################################################################

files = list.files(source.dir,pattern="*.asc")
files2 = list.files(source.dir,pattern="*.gz")
files = setdiff(files,files2)
rm(files2)

setwd(source.dir)

for (file in files) {
  #asc = read.asc.gz(file)
  #ras = raster.from.asc(asc)
  ras = raster(file)
  ext = extent(ras)
  cat("\n\n",ras@file@name,"\n",ext@xmin,ext@xmax,ext@ymin,ext@ymax)
}
