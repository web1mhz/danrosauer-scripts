
library(SDMTools)
source("ClimateByTime.r")

################ set the parameters ################
env.dir         <- "C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim"
model.dir       <- "C:/Users/u3579238/GISData/Helping/Sally/Wyulda_stability/maxent.output_less_bias_2_5deg"
stability.dir   <- "C:/Users/u3579238/GISData/Helping/Sally/Wyulda_stability/stability_2_5deg_less_bias"
stability.file  <- "static.sum.cost.asc"

x               <- 127.82
y               <- -14.61
layer           <- "bioclim_16.asc.gz"
env_value_name  <- "precipitation of wettest qtr"
is_temp         <- FALSE
plot_model      <- TRUE
################ end of parameters ################

# get the climate time series
times <- point.time.series(x,y,layer,env.dir)
cat("\n\n")
times$yearnum <- as.numeric(as.character(times$year))
if (is_temp) {
  times$value <- times$value / 10
}

#get the model time series
if (plot_model) {
  models <- model.time.series(x,y,layer,model.dir)
  models$yearnum <- as.numeric(as.character(models$year))

  #get the stability value
  grid_path     <- paste(stability.dir,stability.file,sep="/")
  stability.asc <- read.asc(grid_path)
  point <- data.frame(x,y)
  stability     <- extract.data(point,stability.asc)
}

#draw the time series graph
windows(10,7)
par(mar=c(5,4,4,5)+.1)

header <- paste("Time series of",env_value_name,"since 120kya\nLatitude:",y,"Longitude:",x)
yrange <- max(times$value) - min(times$value)
ymin <- min(times$value,na.rm=T) - yrange/6
ymax <- max(times$value, na.rm=T) + yrange/6

plot(times$yearnum,times$value,xlim=c(120,0),ylim=c(ymin,ymax),type="l",col="red", lwd=2,main=header,ylab=env_value_name,xlab="thousands of years before present", cex.lab=1.3, cex.axis=1.3)
abline(h=times$value[1],lty=2,col="red",lwd=2)

if (plot_model) {
  model_value_name <- "habitat suitability"
  par(new=TRUE)
  plot(models$yearnum, models$value,,type="l",col="black",xaxt="n",yaxt="n",xlab="",ylab="",lty=2,lwd=2,xlim=c(120,0),ylim=c(0,1.07))
  axis(4,cex.axis=1.3)
  mtext(model_value_name,side=4,line=3, cex=1.5)
  abline(h=stability,lwd=2,col="blue")

  legend("topleft",col=c("red","red","black","blue"),lty=c(1,2,2,1),lwd=2,cex=1.2,legend=c(env_value_name,"current environment",model_value_name,"stability - static"),bg="white")

} else {

  legend("topleft",col=c("red","red"),lty=c(1,2),lwd=2,cex=1.2,legend=c(env_value_name,"current environment"),bg="white")
}

