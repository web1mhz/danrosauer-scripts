rm(list=ls())
library(SDMTools)
source("ClimateByTime.r")

################ set the parameters ################
env.dir         <- "C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim"
model.dir       <- "C:/Users/u3579238/Work/Refugia/Stability/AWT_RF/maxent.output.topo.buf200km"
stability.dir   <- "C:/Users/u3579238/Work/Refugia/Stability/AWT_RF/stability.topo.buf200km"
stability.static.file <- "static.sum.cost.asc"
stability.10m.file    <- "shift.10.asc"
place_name            <- "Lynch's Crater"

x               <- 145.7
y               <- -17.366

layer           <- "bioclim_01.asc.gz"
value_name      <- "AWT rainforest SDM"
env_value_name  <- "Mean Anuual Temperature"
is_temp         <- TRUE
plot_env        <- TRUE
plot_model      <- TRUE

################ end of parameters ################
# get the climate time series
if (plot_env) {
  times <- point.time.series(x,y,layer,env.dir)
  cat("\n\n")
  times$yearnum <- as.numeric(as.character(times$year))
  if (is_temp) {
    times$value <- times$value / 10
  }
}

#get the model time series
if (plot_model) {
  models <- model.time.series(x,y,layer,model.dir)
  models$yearnum <- as.numeric(as.character(models$year))

  #get the stability values
  grid_path     <- paste(stability.dir,stability.static.file,sep="/")
  stability.asc <- read.asc(grid_path)
  point <- data.frame(x,y)
  stability_static     <- extract.data(point,stability.asc)
  grid_path     <- paste(stability.dir,stability.10m.file,sep="/")
  stability.asc <- read.asc(grid_path)
  point <- data.frame(x,y)
  stability_10m     <- extract.data(point,stability.asc)
}

#draw the time series graph
windows(10,8)
par(mar=c(5,4,4,5)+.1)

header <- paste("Time series of",value_name,"since 120kya\nLatitude:",y,"Longitude:",x)
header <- paste(header,"\nModel:",model.dir)

# first plot the model time series
plot(models$yearnum, models$value,,type="l",col="black",xlab="kya",ylab="habitat suitability",lwd=2,xlim=c(120,0),ylim=c(0,1.07), main=header,cex.main=0.9,cex.lab=1.3, cex.axis=1.3)
abline(stability_static,0,lwd=1,lty=3)
abline(stability_10m,0,lwd=1,lty=4)
text(x=60,y=1.075,labels=place_name,cex=1.4)

if (plot_env) {
  yrange <- max(times$value) - min(times$value)
  ymin <- min(times$value,na.rm=T) - yrange/6
  ymax <- max(times$value, na.rm=T) + yrange/6
  
  par(new=TRUE)
  plot(times$yearnum, times$value, type="l", col="green", xaxt="n",yaxt="n",xlab="",ylab="",lwd=1,xlim=c(120,0),ylim=c(ymin,ymax),new=T)
  axis(4,cex.axis=1.3)
  mtext(env_value_name,side=4,line=3, cex=1.5)
  abline(h=times$value[1],lty=2,col="green",lwd=1)

  legend("topleft",col=c("green","green","black","black", "black"),lty=c(1,2,1,3,4),lwd=c(1.5,1,1,1,1),cex=1,legend=c(env_value_name,"current environment",value_name,"stability - static", "stability 10myr"),bg="white")
} else {

  legend("topleft",col="black",lty=c(1,3,4),lwd=c(2,1,1),cex=1.2,legend=c(env_value_name,"stability - static","stability 10m/yr"),bg="white")
}

