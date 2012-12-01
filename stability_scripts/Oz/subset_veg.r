#drafted by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
#GNU General Public License .. feel free to use / distribute ... no warranties

################################################################################
library(SDMTools, raster)

################################################################################
################################################################################

#define directories
work.dir = 'C:/Users/u3579238/Work/Refugia/Stability/AWT_RF/'; setwd(work.dir)
template = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_01.asc'
#veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/Oz/pre1750_mvg.asc.gz'
veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/nvis4_1_mvs_pre_geo_4km_trop_rf_only.asc'

#define some basic data
mvg.asc = read.asc(veg_grid)                 # read in the vegetation grid
rf.asc = mvg.asc; rf.asc[which(is.finite(rf.asc) & rf.asc!=1)] = 0  #set all veg != 1 (rainforests) to 0

################################################################################
#get a subset of the data for occur & background
pos = as.data.frame(which(is.finite(mvg.asc),arr.ind=TRUE)) #get all points that have data
pos$mvg = mvg.asc[cbind(pos$row,pos$col)] #append the vegetation data

mvg.ras = raster.from.asc(mvg.asc)

old_extent = extent(mvg.ras)

ymax=-15.3
ymin=-19.5

smaller_extent = old_extent
smaller_extent@ymax = ymax
smaller_extent@ymin = ymin

template.asc=read.asc(template)
template.ras=raster.from.asc(template.asc)
template.extent=extent(template.ras)

mvg_new.ras = crop(mvg.ras,smaller_extent)
mvg_new.ras = extend(mvg_new.ras,template.extent,0)
mvg_new.asc = asc.from.raster(mvg_new.ras)
mvg_new.asc[which(is.na(template.asc))] = NA


rows_to_include <- y[which(y >= ymin & y <= ymax)]
rownames_to_include <- names(rows_to_include)
rowset <- as.integer(substr(rownames_to_include,2,4))
pos[(pos$row %in% rowset) & (pos$mvg==1),3]


pos.subset = pos[which(pos$row %% 2 == 0 & pos$col %% 2 == 0),] #get a subset of the data as a regular grid of every second cell (even row / col numbers)
#pos.subset= pos  #this line is an alternative to the previous, which selects every 2nd cell
pos.subset = rbind(pos.subset,pos[which(pos$mvg==1),]); pos.subset = unique(pos.subset) #ensure all rf veg cells are included in dataset
