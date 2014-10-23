#drafted by Dan Rosauer

################################################################################
library(SDMTools, raster)

################################################################################

#define directories
work.dir =    'C:/Users/u3579238/GISData/Helping/Rosa/Het_stability'; setwd(work.dir)
template =    'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/landmask.asc'
occ_grid =    'C:/Users/u3579238/GISData/Helping/Rosa/Het_stability/ca6_dist_4deg.asc'
output_grid = 'C:/Users/u3579238/GISData/Helping/Rosa/Het_stability/ca6_dist_4deg_rightextent.asc'
# ymax= -19.3
# ymin= -29
# xmin= 146.5
# xmax= 153.7

ymax= 0
ymin= -90
xmin= 0
xmax= 180

#define some basic data
occ.asc = read.asc(occ_grid)              # read in the occurrence grid
occ.asc[which(is.finite(occ.asc))] <- 1   #set all occurrences to 1

template.asc=read.asc(template)
template.asc[which(is.finite(template.asc))] <- 0

occ.ras = raster.from.asc(occ.asc)
template.ras=raster.from.asc(template.asc)

new_extent = extent(template.ras)
new_extent@ymax = min(ymax,new_extent@ymax)
new_extent@ymin = max(ymin,new_extent@ymin)
new_extent@xmin = max(xmin,new_extent@xmin)
new_extent@xmax = min(xmax,new_extent@xmax)
template_new.ras <- crop(template.ras,new_extent)
occ_new.ras <- crop(occ.ras,new_extent)

occ_out.ras <- merge(occ_new.ras,template_new.ras)
writeRaster(occ_out.ras,output_grid,format="ascii",overwrite=TRUE)
