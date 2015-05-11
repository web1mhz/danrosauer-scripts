rm(list=ls())
library(raster)
library(rgdal)
input_MIROC.dir <- "C:/Users/u3579238/GISData/EnvironmentGrids/Paleo/wc_2_5m_MIROC3.2_21k_bio/2_5m/"
input_CCSM.dir <- "C:/Users/u3579238/GISData/EnvironmentGrids/Paleo/wc_2_5m_CCSM_21k_bio/"
input_HAD3_LGM.dir <- "C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/021"
input_HAD3_now.dir <- "C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000"
output.dir<- "C:/Users/u3579238/GISData/EnvironmentGrids/Paleo/wc_2_5m_MIROC3.2_21k_bio/compare_2_5m/"
output.filename <- paste(output.dir,"CAT2 bioclim CCSM_MIROC_HAD3_NOW.pdf",sep="")
setwd(input_MIROC.dir)

# aus_xmax <- 148
# aus_xmin <- 118.5
# aus_ymax <- -10
# aus_ymin <- -23.5
aus_xmax <- 148
aus_xmin <- 130
aus_ymax <- -11
aus_ymin <- -27

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

file.pattern    <- 'bio+.*asc$'  #regex
setwd(input_HAD3_LGM.dir)
input_files= list.files(path = input_HAD3_LGM.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE,
                        ignore.case = TRUE, include.dirs = FALSE)
lgm_HAD3.stack <- stack(input_files)
lgm_HAD3.stack <- stack(crop(lgm_HAD3.stack,my.extent))
lgm_HAD3.stack <- stack(resample(lgm_HAD3.stack,lgm_ccsm.stack,method="ngb"))
numbers <- c(1,4,10,11,12,15,16,17)
names(lgm_HAD3.stack) <- paste("bio",numbers,sep="")
HAD3_names <- names(lgm_HAD3.stack)

file.pattern    <- 'bio+.*asc$'  #regex
setwd(input_HAD3_now.dir)
input_files= list.files(path = input_HAD3_now.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE,
                        ignore.case = TRUE, include.dirs = FALSE)
now_HAD3.stack <- stack(input_files)
now_HAD3.stack <- stack(crop(now_HAD3.stack,my.extent))
now_HAD3.stack <- stack(resample(now_HAD3.stack,lgm_ccsm.stack,method="ngb"))
numbers <- c(1,4,10,11,12,15,16,17)
names(now_HAD3.stack) <- paste("bio",numbers,sep="")

input.dir <- "C:/Users/u3579238/GISData/EnvironmentGrids/WorldClim/bio_2-5m_bil/"
setwd(input.dir)
file.pattern    <- '*.bil$'  #regex
input_files <- list.files(path = input.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE, ignore.case = TRUE, include.dirs = FALSE)
now.stack <- stack(input_files)
now.stack <- stack(crop(now.stack,my.extent))

#region.shp <- shapefile("C:/Users/u3579238/Work/AMT/Data/Regions/OzMonsoonal/mon_trop_poly_clip.shp")
region.shp <- shapefile("G:/GIS_data/Cat/Cat_region.shp")

# rename lgm layers to match current
names(lgm_miroc.stack) <- names(now.stack)

pdf(output.filename,paper="a4r",width=10.99,height=7.57,onefile=T)

for (name in names(now.stack)) {
  index <- which(names(now.stack)==name)
  HAD3_index <- which(HAD3_names==name)
  include_HAD3 <- name %in% HAD3_names
  
  if (include_HAD3) {

    #windows(18,8)
    par(mfcol=c(3,4))
    
    lgm_MIROC.ras <- lgm_miroc.stack@layers[index][[1]]
    lgm_CCSM.ras  <- lgm_ccsm.stack@layers[index][[1]]
    lgm_HAD3.ras  <- lgm_HAD3.stack@layers[HAD3_index][[1]]
    now_HAD3.ras  <- now_HAD3.stack@layers[HAD3_index][[1]]
    now.ras       <- now.stack@layers[index][[1]]
    
    plot(lgm_MIROC.ras,main=paste("LGM MIROC",name))
    plot(region.shp,add=T)
    
    plot(lgm_CCSM.ras,main=paste("LGM CCSM",name))  
    plot(region.shp,add=T)
    
    plot(lgm_HAD3.ras,main=paste("LGM HAD3",name))  
    plot(region.shp,add=T)
    
    plot(now.ras,main=paste("Now",name))
    plot(region.shp,add=T)
    
    diff.ras <- now.ras - lgm_MIROC.ras
    plot(diff.ras,main=paste("Diff (now - MIROC LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)
    
    diff.ras <- now.ras - lgm_CCSM.ras
    plot(diff.ras,main=paste("Diff (now - CCSM LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)
    
    diff.ras <- now.ras - lgm_HAD3.ras
    plot(diff.ras,main=paste("Diff (now - HAD3 LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)

    diff.ras <- now_HAD3.ras - lgm_HAD3.ras
    plot(diff.ras,main=paste("Diff (HAD3 now - HAD3 LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)
        
    diff_model.ras <- lgm_MIROC.ras - lgm_CCSM.ras
    plot(diff_model.ras,main=paste("Diff (MIROC LGM - CCSM LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)
    
    diff_model.ras <- lgm_HAD3.ras - lgm_CCSM.ras
    plot(diff_model.ras,main=paste("Diff (HAD3 LGM - CCSM LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)
    
    diff_model.ras <- lgm_HAD3.ras - lgm_MIROC.ras
    plot(diff_model.ras,main=paste("Diff (HAD3 LGM - MIROC LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)

  } else {
  
    par(mfcol=c(2,3))
      
    lgm_MIROC.ras <- lgm_miroc.stack@layers[index][[1]]
    lgm_CCSM.ras  <- lgm_ccsm.stack@layers[index][[1]]
    if (include_HAD3) {
      lgm_HAD3.ras <- lgm_HAD3.stack@layers[index][[1]]
    }
    now.ras <- now.stack@layers[index][[1]]
    
    plot(lgm_MIROC.ras,main=paste("LGM MIROC",name))
    plot(region.shp,add=T)
    
    plot(lgm_CCSM.ras,main=paste("LGM CCSM",name))  
    plot(region.shp,add=T)
    
    plot(now.ras,main=paste("Now",name))
    plot(region.shp,add=T)
    
    diff.ras <- now.ras - lgm_MIROC.ras
    plot(diff.ras,main=paste("Diff (now - MIROC LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)
    
    diff.ras <- now.ras - lgm_CCSM.ras
    plot(diff.ras,main=paste("Diff (now - CCSM LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)
    
    diff_model.ras <- lgm_MIROC.ras - lgm_CCSM.ras
    plot(diff_model.ras,main=paste("Diff (MIROC LGM - CCSM LGM)",name),col=heat.colors(255))
    plot(region.shp,add=T)
  }
}
dev.off()
