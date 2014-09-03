# this script will extract records from the Moritz Lab database and generate maps and kml files
rm(list=ls())

library(maptools)
library(sp)
library(raster)

# parameters
output_prefix <- 'AMT narrow2'
models_dir <- "C:/Users/u3579238/Work/Refugia/Stability/Kimberley/maxent.output_narrow_AMT2/"
output_dir <- "C:/Users/u3579238/Work/Refugia/Stability/Kimberley/maps_narrow_AMT2"
per_pic   <- 4

subset <- c(1, 4, 7, 13, 19, 23, 27, 32, 42, 52, 57, 62)
do_subset=T
#################

setwd(models_dir)

models <- dir(,path=models_dir,pattern = ".asc",include.dirs = T,full.names = T)
if (do_subset) {models <- models[subset]}

coast.shp <- shapefile("C:/Users/u3579238/GISData/aus_1m.shp")
coast.shp <- coast.shp[coast.shp$NAM != "PAPUA NEW GUINEA" & coast.shp$NAM != "INDONESIA",]

IBRA.shp <- shapefile("C:/Users/u3579238/GISData/IBRA/IBRA7_regions.shp")
#region_list <- "Northern Kimberley"
#region_list <- c("Northern Kimberley", "Central Kimberley", "Victoria Bonaparte")
#region_list <- c("Arnhem Coast","Arnhem Plateau","Central Arnhem","Daly Basin","Darwin Coastal","Gulf Fall and Uplands",
#                  "Gulf Coastal","Mount Isa Inlier","Ord Victoria Plain","Pine Creek","Sturt Plateau","Tiwi Cobourg",
#                  "Northern Kimberley", "Central Kimberley", "Victoria Bonaparte")
region_list <- c("Arnhem Coast","Arnhem Plateau","Central Arnhem","Daly Basin","Darwin Coastal","Gulf Fall and Uplands",
                 "Gulf Coastal","Mount Isa Inlier","Pine Creek","Tiwi Cobourg","Northern Kimberley", "Central Kimberley",
                 "Victoria Bonaparte")
IBRA.shp    <- IBRA.shp[IBRA.shp$REG_NAME_7 %in% region_list,]


pic_count <- 0

for (model in models) {

  cat("mapping ",model,"\n")
  model.ras  <- raster(model)
  model_date <- substr(model,nchar(model)-6,nchar(model)-4)
  if (do_subset) {
    filename <- paste(output_dir,"/",output_prefix," selected_",model_date,"_x",per_pic,".png",sep="")
  } else {
    filename <- paste(output_dir,"/",output_prefix,"_",model_date,"_x",per_pic,".png",sep="")
  }
  # make a map
  #windows(19.2,12)

  if (pic_count %% per_pic ==0) {
    png(filename,width=1920,height=1200)
    par(mfrow=c(2,2))
  }

  main_header <- paste(output_prefix ,model_date, sep=" ")

  plot(model.ras, main = main_header, xlab="longitude", ylab="latitude", ylim=c(-24,-10), xlim = c(118,142), cex.main=1.5)
  plot(coast.shp,lwd=0.6, add=T)
  plot(IBRA.shp,lwd=1.5, add=T)

  if (pic_count %% per_pic == per_pic - 1) {
    dev.off()
  }
  pic_count <- pic_count + 1
}

graphics.off()
