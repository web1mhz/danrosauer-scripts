# a script to summarize values from one ascii grid to a coarser (but not necessarily aligned) grid
# Dan Rosauer - May 2013

library(raster)

# functions to create a circular filter for
make_circ_filter<-function(radius, res){
  circ_filter<-matrix(NA, nrow=1+(2*radius/res), ncol=1+(2*radius/res))
  dimnames(circ_filter)[[1]]<-seq(-radius, radius, by=res)
  dimnames(circ_filter)[[2]]<-seq(-radius, radius, by=res)
  sweeper<-function(mat){
    for(row in 1:nrow(mat)){
      for(col in 1:ncol(mat)){
        dist<-sqrt((as.numeric(dimnames(mat)[[1]])[row])^2 +
                     (as.numeric(dimnames(mat)[[1]])[col])^2)
        if(dist<=radius) {mat[row, col]<-1}
      }
    }
    return(mat)
  }
  out<-sweeper(circ_filter)
  return(out)
}

########## Parameters ##########
input.raster      <- 'C:/Users/u3579238/GISData/Mackey/Greenspots_Apr13/mean00_12.asc'
template.raster   <- 'C:/Users/u3579238/GISData/EnvironmentGrids/AusGDMGrids/ForMaxent/bio1.asc'
output.raster     <- 'C:/Users/u3579238/GISData/Mackey/Greenspots_Apr13/coarser/mean00_12_max1.5km_01deg.asc'
funct             <- 'max'
radius.cells      <- 5
################################

input.ras           <- raster(input.raster)
template.ras        <- raster(template.raster)
projection(template.ras) <- projection(input.ras)

template.ext        <- extent(template.ras)
template.ext@xmin   <- 120.5
template.ext@xmax   <- 132
template.ext@ymax   <- -12.5
template.ext@ymin   <- -20.5

crop.template.ras   <- crop(template.ras,template.ext)
crop.input.ras      <- crop(input.ras,template.ext)

resol               <- xres(crop.input.ras)
filter              <- make_circ_filter(radius=(resol*radius.cells),resol)

focal.ras     <- focal(crop.input.ras,w=filter,fun=max,na.rm=TRUE)
focal.resample.ras      <- resample(focal.ras,crop.template.ras,method="bilinear")

windows(8,8)
plot(crop.input.ras)
windows(8,8)
plot(focal.resample.ras)

writeRaster(focal.resample.ras,output.raster,format="ascii",overwrite=TRUE)
