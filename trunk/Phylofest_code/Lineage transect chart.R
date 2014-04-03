rm(list=ls())

library(raster)
library(sp)
library(stringr)
library(ggplot2)
library(reshape2)

coords <- function(the.line,n,position){
  line_length <- LineLength(the.line)
  dist_from_start <- position * (line_length / n)
  x <- the.line@coords[1,1] + (the.line@coords[2,1] - the.line@coords[1,1]) * (position/n)
  y <- the.line@coords[1,2] + (the.line@coords[2,2] - the.line@coords[1,2]) * (position/n)
  return(as.data.frame(cbind(x,y)))
}


#define directories
input.dir       <- 'C:/Users/u3579238/Work/AMT/Models/lineage_models/asc_clipped_cube_method/'
output.dir      <- 'C:/Users/u3579238/Work/AMT/Maps/lineage_models/Heteronotia/'
file_pattern    <- '[:alnum:]*Heteronotia_b.*asc$'
group_lin_file  <- 'C:/Users/u3579238/Work/AMT/Models/group_lineage_list.csv'
preface         <- 'lin_model_'
suffix          <- ''

lineage_site_file   <- 'C:/Users/u3579238/Work/AMT/Models/lineage_sites/Heteronotia_lin_loc_from_db_18mar14_nth_of_22S.csv'

too.low = 0

line_coords <- data.frame(cbind(long=c(133,134),lat=c(-11.5,-22)))
#line_coords <- data.frame(cbind(long=c(127.4,118.8),lat=c(-14,-21.8)))
#line_coords <- data.frame(cbind(long=c(121.8,132.5),lat=c(-19,-19)))

####  end of parameters

files = list.files(path = input.dir, pattern = file_pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

setwd(input.dir)

lin.stack <- stack(files)

i <- 1
lineage <- vector(mode="character")

splitfiles <- strsplit(files,"_")
for (file in splitfiles) {
  this.lineage <- file[5]
  lineage[i] <- strsplit(this.lineage,"[.]")[[1]][1]
  i <- i+1
}

names(lin.stack) <- lineage

# make the spatial lines object
the.line <- Line(line_coords)
the.lines<- Lines(list(the.line),"myLine")
the.spatiallines <- SpatialLines(list(the.lines))
rm(the.lines)

# make a spatial lines data frame
data <- as.data.frame(1,row.names="myLine")
the.spatiallinesdf <- SpatialLinesDataFrame(the.spatiallines,data)

#write the line as a shapefile
#shapefile(the.spatiallinesdf,paste(output.dir,"line3.shp",sep=""),overwrite=T)

lin_line_values <- extract(lin.stack,the.spatiallines)
lin_line_values <- lin_line_values[[1]]

#group low scoring lineages as 'other'
lin_max <- apply(lin_line_values,MARGIN=2,FUN=max,na.rm=T)
minors <- which(lin_max < too.low)
main_values <- lin_line_values[,-minors]
n <- ncol(main_values)
other <- apply(lin_line_values[,minors],MARGIN=1,FUN=sum)

rows <- nrow(lin_line_values)
the_coords <- round(coords(the.line,rows,1:rows),3)
row.names(main_values) <- the_coords[,1]
#y_coords <- the_coords[main_values$position,2]

main_values <- cbind(main_values,other)
main_values <- melt(main_values,varnames=c("position","lineage"))

windows()
ggplot(main_values, aes(x=position, y=value, fill=lineage)) + geom_area() + ggtitle("Lineage values along a transect")


