#drafted by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
#GNU General Public License .. feel free to use / distribute ... no warranties

################################################################################
library(SDMTools, raster)

################################################################################
################################################################################

#define directories
template = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/meq_rainforest2.asc'
#veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/Oz/pre1750_mvg.asc.gz'
veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/nvis4_1_mvs_pre_geo_4km_trop_rf_only.asc'
output_grid = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/cec_rainforest.asc'

#define some basic data
mvg.asc = read.asc(veg_grid)                 # read in the vegetation grid
rf.asc = mvg.asc; rf.asc[which(is.finite(rf.asc) & rf.asc!=1)] = 0  #set all veg != 1 (rainforests) to 0

################################################################################
#get a subset of the data for occur & background
pos = as.data.frame(which(is.finite(mvg.asc),arr.ind=TRUE)) #get all points that have data
pos$mvg = mvg.asc[cbind(pos$row,pos$col)] #append the vegetation data

mvg.ras = raster.from.asc(mvg.asc)

old_extent = extent(mvg.ras)

ymax=-23.8
ymin=-33.0

smaller_extent = old_extent
smaller_extent@ymax = ymax
smaller_extent@ymin = ymin

template.asc=read.asc(template)
template.ras=raster.from.asc(template.asc)
template.extent=extent(template.ras)

mvg_new.ras = crop(mvg.ras,smaller_extent)
mvg_new.ras = extend(mvg_new.ras,template.extent,0)
mvg_new.ras = crop(mvg_new.ras,template.extent)
mvg_new.asc = asc.from.raster(mvg_new.ras)
mvg_new.asc[which(is.na(template.asc))] = NA

write.asc(mvg_new.asc,output_grid)
