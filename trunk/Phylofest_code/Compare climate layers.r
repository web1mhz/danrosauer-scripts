rm(list=ls())
library(raster)
library(rgdal)
input_MIROC.dir <- "C:/Users/u3579238/GISData/EnvironmentGrids/Paleo/wc_2_5m_MIROC3.2_21k_bio/2_5m/"
input_CCSM.dir <- "C:/Users/u3579238/GISData/EnvironmentGrids/Paleo/wc_2_5m_CCSM_21k_bio/"
output.dir<- "C:/Users/u3579238/GISData/EnvironmentGrids/Paleo/wc_2_5m_MIROC3.2_21k_bio/compare_2_5m/"
setwd(input_MIROC.dir)

aus_xmax <- 148
aus_xmin <- 118.5
aus_ymax <- -10
aus_ymin <- -24

my.extent <- extent(c(aus_xmin,aus_xmax,aus_ymin,aus_ymax))
rm(aus_xmax,aus_xmin,aus_ymax,aus_ymin)

file.pattern    <- '*.bil$'  #regex
input_files= list.files(path = input_MIROC.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE,
                        ignore.case = TRUE, include.dirs = FALSE)
lgm_miroc.stack <- stack(input_files)
lgm_miroc.stack <- stack(crop(lgm_miroc.stack,my.extent))

setwd(input_CCSM.dir)
input_files= list.files(path = input_CCSM.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE,
                        ignore.case = TRUE, include.dirs = FALSE)
lgm_ccsm.stack <- stack(input_files)
lgm_ccsm.stack <- stack(crop(lgm_ccsm.stack,my.extent))

input.dir <- "C:/Users/u3579238/GISData/EnvironmentGrids/WorldClim/bio_2-5m_bil/"
setwd(input.dir)
input_files <- list.files(path = input.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE, ignore.case = TRUE, include.dirs = FALSE)
now.stack <- stack(input_files)
now.stack <- stack(crop(now.stack,my.extent))

region.shp <- shapefile("C:/Users/u3579238/Work/AMT/Data/Regions/OzMonsoonal/mon_trop_poly_clip.shp")

# rename lgm layers to match current
names(lgm_miroc.stack) <- names(now.stack)

for (name in names(now.stack)) {
  filename <- paste(output.dir,"CCSM_MIROC_NOW_",name,".pdf",sep="")
  #png(filename,width=18,height=8,units="in",res=150)
  pdf(filename,width=18,height=8)
  #windows(18,6)
  par(mfcol=c(2,3))
  
  index <- which(names(now.stack)==name)
  
  lgm_MIROC.ras <- lgm_miroc.stack@layers[index][[1]]
  lgm_CCSM.ras  <- lgm_ccsm.stack@layers[index][[1]]  
  now.ras <- now.stack@layers[index][[1]]
  
  plot(lgm_MIROC.ras,main=paste("LGM MIROC",name))
  plot(region.shp,add=T)
  
  plot(lgm_CCSM.ras,main=paste("LGM CCSM",name))  
  plot(region.shp,add=T)
  
  plot(now.ras,main=paste("Now",name))
  plot(region.shp,add=T)
  
  diff.ras <- now.ras - lgm_MIROC.ras
  plot(diff.ras,main=paste("Difference (now - MIROC LGM)",name))
  plot(region.shp,add=T)
  
  diff.ras <- now.ras - lgm_CCSM.ras
  plot(diff.ras,main=paste("Difference (now - CCSM LGM)",name))
  plot(region.shp,add=T)
  
  diff_model.ras <- lgm_MIROC.ras - lgm_CCSM.ras
  plot(diff_model.ras,main=paste("Difference (MIROC LGM - CCSM LGM)",name))
  plot(region.shp,add=T)
  dev.off() 
}
