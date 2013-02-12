#drafted by Dan Rosauer, using elements from a script by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )

################################################################################
library(SDMTools, raster)

################################################################################
################################################################################

#define directories
template = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/landmask.asc'
#veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/Oz/pre1750_mvg.asc.gz'
veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/nvis4_1_mvs_pre_geo_4km_rf_only.asc'
output_grid = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/CEC_all_rainforest.asc'

#define some basic data
mvs.asc = read.asc(veg_grid)                 # read in the vegetation grid
template.asc=read.asc(template)
rf.asc = mvs.asc 
#rf.asc[which(is.finite(rf.asc) & !(rf.asc==1 | rf.asc==2 | rf.asc==6))] <- 0  #set all veg != 1 (rainforests) to 0
rf.asc[which(is.finite(rf.asc))] <- 1  #set all rainforest to 1

################################################################################
#get a subset of the data for occur & background
pos = as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data
pos$rf = rf.asc[cbind(pos$row,pos$col)] #append the vegetation data

rf.ras = raster.from.asc(rf.asc)

old_extent = extent(rf.ras)

ymax=-23.8
ymin=-32.5

smaller_extent = old_extent
smaller_extent@ymax = ymax
smaller_extent@ymin = ymin

template.ras=raster.from.asc(template.asc)
template.extent=extent(template.ras)

rf_new.ras = crop(rf.ras,smaller_extent)
rf_new.ras = extend(rf_new.ras,template.extent,0)
rf_new.ras = crop(rf_new.ras,template.extent)
rf_new.asc = asc.from.raster(rf_new.ras)
rf_new.asc[which(is.na(template.asc))] = NA

rf_new.asc[is.na(rf_new.asc) & !is.na(template.asc)] <- 0

write.asc(rf_new.asc,output_grid)
cat("\nWritten",output_grid,"\n")
