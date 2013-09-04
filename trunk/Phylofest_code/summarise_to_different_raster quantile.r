# a script to summarize values from one ascii grid to a coarser grid
# Dan Rosauer - August 2013

rm(list=ls())
library(raster)

my.quantile <- function(values,quant) {
  out <- as.numeric(quantile(values,quant,na.rm=T))
  return(out)
}

########## Parameters ##########
input.raster      <- 'C:/Users/u3579238/GISData/Helping/Sally/Wyulda_models/slope.asc'
template.raster   <- 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_01.asc'
output.raster     <- 'C:/Users/u3579238/GISData/Helping/Sally/median.asc'
################################

input.ras           <- raster(input.raster)
template.ras        <- raster(template.raster)
projection(template.ras) <- projection(input.ras)

template.ext        <- extent(template.ras)
# template.ext@xmin   <- 120.5
# template.ext@xmax   <- 132
# template.ext@ymax   <- -12.3
# template.ext@ymin   <- -20.8

# crop rasters to the same extent
crop.template.ras   <- crop(template.ras,template.ext)
ext.input.ras       <- extend(input.ras,template.ext)
crop.input.ras      <- crop(input.ras,template.ext)

# create a zonal raster based on the cells in the template raster
non.NA.length <- length(which(!is.na(crop.template.ras[])))
crop.template.ras[which(!is.na(crop.template.ras[]))] <- rep(1:non.NA.length)
zones.ras   <-        resample(x=crop.template.ras, crop.input.ras, method="ngb")

#zonal.result     <- zonal(crop.input.ras,zones.ras,digits=8,fun="my.quantile",quant=0.75)
zonal.result     <- as.data.frame(zonal(crop.input.ras,zones.ras,digits=8,fun=median,na.rm=T))
names(zonal.result) <- c("zone","value")
zones  <- data.frame(crop.template.ras[])
zones <- cbind(1:nrow(zones),zones)
names(zones) <- c("cellID","zone")
result <- merge(zones,zonal.result,by="zone",all.x=T,all.y=F)
result  <- result[order(result$cellID),]

result.ras       <- crop.template.ras
result.ras[]      <- result[,3]

windows(8,8)
plot(crop.input.ras, main=paste("Input:\n",input.raster))
windows(8,8)
plot(result.ras, main=paste("Output:\n",output.raster))

result.ras
writeRaster(result.ras,output.raster,format="ascii",overwrite=TRUE)

#now redo this for elevation range
########## Parameters ##########
input.raster      <- 'C:/Users/u3579238/GISData/EnvironmentGrids/Topo/etopo1_landsea/etopo1_aust_and_seas.asc'
output.raster     <- 'C:/Users/u3579238/GISData/Helping/Sally/elev_range_bathym.asc'
################################

#input.ras           <- raster(input.raster)

# The following lines, replace the raster file input above, as SDM tools is more forgiving of minor errors
library(SDMTools)
input.asc <- read.asc(input.raster)
input.ras <- raster.from.asc(input.asc)
rm(input.asc)
# comment out the section above for normal, working rasters

my.crs <- CRS("+proj=longlat +datum=WGS84")
projection(input.ras) <- my.crs

# crop rasters to the same extent
ext.input.ras       <- extend(input.ras,template.ext)
crop.input.ras      <- crop(input.ras,template.ext)

zones.ras   <-        resample(x=crop.template.ras, crop.input.ras, method="ngb")

zonal.result.max     <- as.data.frame(zonal(crop.input.ras,zones.ras,digits=8,fun=max,na.rm=T))
zonal.result.max[which(abs(zonal.result.max[2])==Inf),2]  <- NA
zonal.result.min     <- as.data.frame(zonal(crop.input.ras,zones.ras,digits=8,fun=min,na.rm=T))
zonal.result.min[which(abs(zonal.result.min[2])==Inf),2]  <- NA
zonal.result   <- as.data.frame(cbind(zonal.result.max[,1],zonal.result.max[,2] - zonal.result.min[,2]))
#rm(zonal.result.max,zonal.result.min)

names(zonal.result) <- c("zone","value")
zones  <- data.frame(crop.template.ras[])
zones <-  cbind(1:nrow(zones),zones)
names(zones) <- c("cellID","zone")
result <- merge(zones,zonal.result,by="zone",all.x=T,all.y=F)
result  <- result[order(result$cellID),]

result.ras       <- crop.template.ras
result.ras[]      <- result$value

windows(8,8)
plot(crop.input.ras, main=paste("Input:\n",input.raster))
windows(8,8)
plot(result.ras, main=paste("Output:\n",output.raster))

result.ras
writeRaster(result.ras,output.raster,format="ascii",overwrite=TRUE)

###################################
# calculate slope for this raster # 
output.raster     <- 'C:/Users/u3579238/GISData/Helping/Sally/slope_median_bathym.asc'

slope_bathym.ras <- terrain(crop.input.ras,opt="slope",unit="degrees")

zonal.result     <- as.data.frame(zonal(slope_bathym.ras,zones.ras,digits=8,fun=median,na.rm=T))
names(zonal.result) <- c("zone","value")
zones  <- data.frame(crop.template.ras[])
zones <- cbind(1:nrow(zones),zones)
names(zones) <- c("cellID","zone")
result <- merge(zones,zonal.result,by="zone",all.x=T,all.y=F)
result  <- result[order(result$cellID),]

result.ras       <- crop.template.ras
result.ras[]      <- result[,3]

windows(8,8)
plot(crop.input.ras, main=paste("Input:\n",input.raster))
windows(8,8)
plot(result.ras, main=paste("Output:\n",output.raster))

result.ras
writeRaster(result.ras,output.raster,format="ascii",overwrite=TRUE)

