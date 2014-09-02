#drafted by Dan Rosauer, using elements from a script by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )

################################################################################
rm(list=ls())
library(SDMTools)
library(raster)

#define directories
template = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_01.asc'
#veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/Oz/pre1750_mvg.asc.gz'
veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/nvis4_1_mvs_pre_geo_4km_rf_only.asc'
output_suffix = '_all_rainforest_100km'
output_path = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/'
buffer_width = 100000 # background region around RF pixels, in metres

region_lat_boundaries <- c(-9.75,-14.85,-19.6,-23.4,-32.44, -44.3)
region_names          <- c("CYP","AWT","MEQ","CEC","SEA")

#################################################################
template.ras=raster(template)
template.extent=extent(template.ras)

for (region in 1:length(region_names)) {
  ymax = region_lat_boundaries[region]
  ymin = region_lat_boundaries[region+1]
  region_name = region_names[region]

  output_veg  <- paste(output_path,region_name,output_suffix,'.asc',sep='')
  output_mask <- paste(output_path,region_name,output_suffix,'_mask.asc',sep='')

  
  mvs.asc = read.asc(veg_grid)                 # read in the vegetation grid
  rf.asc = mvs.asc 
  #rf.asc[which(is.finite(rf.asc) & !(rf.asc==1 | rf.asc==2 | rf.asc==6))] <- 0  #set all veg != 1 (rainforests) to 0
  rf.asc[which(is.finite(rf.asc))] <- 1  #set all rainforest to 1

  #get a subset of the data for occur & background
  pos = as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data
  pos$rf = rf.asc[cbind(pos$row,pos$col)] #append the vegetation data
  rf.ras = raster.from.asc(rf.asc)
  
  old_extent = extent(rf.ras)
  
  smaller_extent = old_extent
  smaller_extent@ymax = ymax
  smaller_extent@ymin = ymin
  
  rf_new.ras = crop(rf.ras,smaller_extent)
  rf_new.ras = extend(rf_new.ras,template.extent,0)
  rf_new.ras = crop(rf_new.ras,template.extent)
  rf_new.ras[is.na(template.ras)] = NA
  rf_new.ras[is.na(rf_new.ras) & !is.na(template.ras)] <- 0
  
  writeRaster(rf_new.ras,output_veg,overwrite=TRUE)
  cat("\nWritten",output_veg,"\n")
  
  rm(rf.asc, rf.ras)
  
  # now create a buffered mask for this rf region
  region.ras <- raster(output_veg)
  region.ras[region.ras==0] <- NA
  region_buf.ras <- buffer(region.ras,width=buffer_width, doEdge=TRUE)
  region_buf.ras <- mask(region_buf.ras,template.ras)
  writeRaster(region_buf.ras,output_mask,overwrite=TRUE)
  cat("Written",output_mask,"\n")
  
  windows(8,8)
  plot(region_buf.ras,xlim=c(140,160))
  
  Sys.sleep(1)
  
}