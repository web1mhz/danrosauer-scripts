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
input.raster      <- 'C:/Users/u3579238/GISData/Mackey/Greenspots_Apr13/fgreen00_12_yr_min.flt'
#input.raster      <- 'C:/Users/u3579238/GISData/EnvironmentGrids/FractionalCover/fractcover.v2_2.2012.arch.mean.bs_01.asc'
template.raster   <- 'F:/Transfer/env_rasters/adefi_amt_rept_2deg_msk.flt'
output.raster     <- 'C:/Users/u3579238/GISData/Mackey/Greenspots_Apr13/fgreen00_12_yr_min_max005deg_AMT.flt'
funct             <- 'max'
radius.degrees      <- 0.005
################################

input.ras           <- raster(input.raster)
template.ras        <- raster(template.raster)
projection(template.ras) <- projection(input.ras)

template.ext        <- extent(template.ras)
# template.ext@xmin   <- 120.5
# template.ext@xmax   <- 132
# template.ext@ymax   <- -12.5
# template.ext@ymin   <- -20.5
#crop.template.ras   <- crop(template.ras,template.ext)

crop.input.ras      <- crop(input.ras,template.ext)

resol               <- xres(crop.input.ras)
filter              <- make_circ_filter(radius=radius.degrees,resol)

focal.ras     <- focal(crop.input.ras,w=filter,fun=max,na.rm=TRUE)
focal.resample.ras      <- resample(focal.ras,crop.template.ras,method="bilinear")

template_null <- which(is.na(template.ras[]))
focal.resample.ras[template_null] <- NA

windows(8,8)
plot(crop.input.ras)
windows(8,8)
plot(focal.resample.ras)

writeRaster(focal.resample.ras,output.raster,format="ascii",overwrite=TRUE)

# now check that all cells in the template have valid data in the new layer
template_valid <- which(!is.na(template.ras[]))
new_null    <- which(is.na(focal.resample.ras[]))
new_gaps       <- intersect(template_valid,new_null)
new_gaps_orig <- new_gaps

resol  <- res(template.ras)[1]
radius <- 2
gap_filter              <- make_circ_filter(radius=(resol * radius),resol)

i <- 1

orig.resample.ras <- focal.resample.ras
previous_gaps <- length(new_gaps)

while (length(new_gaps) > 0) { # repeat to fill remaining holes
  new_smoothed.ras <- focal.resample.ras

  cat("smoothing run number:", i, "\tgap pixels:", length(new_gaps), "\tradius:", radius, "\n")
  new_smoothed.ras <- focal(focal.resample.ras,w=gap_filter,fun=median,na.rm=TRUE)
  focal.resample.ras[new_gaps] <- new_smoothed.ras[new_gaps]

  new_null    <- which(is.na(new_smoothed.ras[]))
  new_gaps       <- intersect(template_valid,new_null)

  focal.resample.ras[new_gaps] <- new_smoothed.ras[new_gaps_orig]

  if (length(new_gaps) == previous_gaps) {
    radius <- radius + 1
    gap_filter <- make_circ_filter(radius=(resol * radius),resol)
    cat("Expanding radius to:",radius,"\n")
  }
  previous_gaps <- length(new_gaps)
  i <- i + 1
}

new_final.ras <- focal.resample.ras
new_final.ras[new_gaps_orig] <- new_smoothed.ras[new_gaps_orig]
writeRaster(focal.resample.ras,output.raster,format="ascii",overwrite=TRUE)
