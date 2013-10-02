# map Global Mammal Phylogenetic output files for Marxan
# Dan Rosauer - started October 2013

rm(list=ls())
library(gdata)

base.dir            <- "C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/"
input.dir           <- "C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/input/"
output.dir          <- "C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/output2/"

pu_list.file        <- paste(input.dir,"pu_id_lookup.csv",sep="")
marxan_result.file  <- paste(output.dir,"output_ssoln.txt",sep="")

# mapping options
draw_map <- TRUE
draw_raster <- TRUE  #if false, plot points

############################################
setwd(base.dir)

pu_list <- read.csv(pu_list.file)
marxan_result <- read.csv(marxan_result.file)
max_count <- max(marxan_result$number)

marxan_xy <- merge(pu_list,marxan_result,by.x="pu_id",by.y="planning_unit")
marxan_xy$X <- NULL

setwd(output.dir)
write.csv(marxan_xy,"marxan_xy.csv")


if (draw_map) {
  # make a map
  library(classInt)
  library(maptools)
  
  class_count <- 10
  my.classes  <- classIntervals(marxan_xy$number,n=class_count,style="equal",digits=2)
  my.pal<-c("lightgrey","yellow","red","purple") # choose colors
  
  windows()
  header="Selecting areas which maximise representation of mammal PD"
  
  if (draw_raster) {
    library(raster)
    my.colorRampPalette <- colorRampPalette(my.pal)
    marxan.ras <- rasterFromXYZ(marxan_xy[,3:5])
    plot(marxan.ras,col=my.colorRampPalette(class_count),ylim=c(0,140),main=header)
  } else {
    my.col <- findColours(my.classes,my.pal,over=">",under="<",cutlabels=F)
    legtext <- names(attr(my.col,"table"))  # declare labels
    legcols <- attr(my.col,"palette")
    
    plot(marxan_xy$Axis_0,marxan_xy$Axis_1,col=my.col,pch=20, main=header)
    legend(0,60,legend=legtext,col=my.col,fill=legcols)
  }

  #world <- readShapeLines("C:/Users/u3579238/GISData/ne_50m_coastline/ne_50m_coastline.shp")
  #plot(world, lwd=0.5, add=T, col="black")

}