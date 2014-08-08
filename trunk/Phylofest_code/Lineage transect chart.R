

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

moving_avg <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

#######################################################
##############  Main Program Starts Here  #############
#######################################################

#define directories
input.dir       <- 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models_aligned_rf_strict/'
output.dir      <- 'C:/Users/u3579238/Work/Phylofest/Maps/'
file_pattern    <- '*Saproscincus_r.*asc$'
#group_lin_file  <- 'C:/Users/u3579238/Work/AMT/Models/group_lineage_list.csv'

smooth_width = 12  # use 0 for no smoothing
too.low = 0.1

#line_coords <- data.frame(cbind(long=c(133,134),lat=c(-11.5,-22)))
#line_coords <- data.frame(cbind(long=c(127.4,118.8),lat=c(-14,-21.8)))
#line_coords <- data.frame(cbind(long=c(121.8,132.5),lat=c(-19,-19)))
#line_coords <- data.frame(cbind(long=c(134,130),lat=c(-12.3,-18)))
line_coords <- data.frame(cbind(long=c(152.92,151.75),lat=c(-26.47,-32.5)))

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
shapefile(the.spatiallinesdf,paste(output.dir,"rosei_line.shp",sep=""),overwrite=T)

lin_line_values <- extract(lin.stack,the.spatiallines)
lin_line_values <- lin_line_values[[1]]

#group low scoring lineages as 'other'
lin_max <- apply(lin_line_values,MARGIN=2,FUN=max,na.rm=T)
minors <- which(lin_max < too.low)
if (length(minors) > 0) {
  main_values <- lin_line_values[,-minors]
  other <- apply(lin_line_values[,minors],MARGIN=1,FUN=sum)
  main_values <- cbind(main_values,other)
} else {
  main_values <- lin_line_values[,]
}

n <- ncol(main_values)

# smooth the values
if (smooth_width>0) {
  main_values_orig <- main_values
  main_values_smooth <- main_values
  for (i in 1:n) {
    main_values_smooth[,i] <- moving_avg(main_values[,i],smooth_width, centered = TRUE)
  }
  main_values <- main_values_smooth
  rm(main_values_smooth)
}

rows <- nrow(lin_line_values)
the_coords <- round(coords(the.line,rows,1:rows),3)
row.names(main_values) <- the_coords[,2]  # column 2, to use latitude as the x axis on the plot in this case
#y_coords <- the_coords[main_values$position,2]

main_values <- melt(main_values,varnames=c("position","lineage"))

title <- "Lineage values along a transect"
if (smooth_width > 0) {
  title <- paste(title,"\nsmooth width =",smooth_width)
}

main_values <- main_values[order(main_values$position,decreasing = F),]

windows(6,12)

my.plot <- ggplot(main_values, aes(x=position, y=value, fill=lineage)) + geom_area() + ggtitle(title) + xlim(min(main_values$position),max(main_values$position)) + coord_flip()
my.plot <- my.plot + labs(y = "Predicted suitability", x = "Latitude") + theme_bw(base_size = 12) + theme(panel.grid.major = element_blank()) + theme(panel.border = element_blank())
print(my.plot)



