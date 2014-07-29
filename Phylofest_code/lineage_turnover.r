rm(list=ls())

library(ecodist)
library(MASS)

setwd("C:/Users/u3579238/Work/Refugia/Results/outputs_jul2014")
models <- read.csv("rf_lin_models.csv")
resolution <- 0.2

border <- "C:/Users/u3579238/GISData/aus5m_coast_no_norfolk.shp"

#strip out empty cells
models <- models[which(models$total > 0),]

#group cells
models$lat_group <- as.integer(models$lat/resolution) * resolution
models$long_group <- as.integer(models$long/resolution) * resolution
models_only <- cbind(models$lat_group,models$long_group,paste(models$lat_group,"_",models$long_group,sep=""),models[,12:(ncol(models)-4)])
names(models_only)[1:3] <- c("lat_group","long_group","cell_name")

models_grouped <- aggregate(models_only[,4:ncol(models_only)],by=list(models_only$lat_group,models_only$long_group,models_only$cell_name),FUN=sum)
names(models_grouped)[1:3] <- c("lat_group","long_group","cell_name")

# set NA values to 0
for (col in 1:ncol(models_grouped)) {
  na_count <- length(which(is.na(models_grouped[,col])))
  cat(names(models_grouped)[col], "\t",na_count,"\n")
  if (na_count>0) {
    models_grouped[which(is.na(models_grouped[,col])),col] <- 0
  }
}

row.names(models_grouped) <- models_grouped$cell_name

dist.dist <- distance(as.matrix(models_grouped[,4:ncol(models_grouped)]), method="bray-curtis",spweight="absence")
dist.dist[which(is.na(dist.dist))] <- 1
dist_matrix <- as.matrix(dist.dist)
gc()

clust <- hclust(dist.dist)

dist_matrix[which(is.na(dist_matrix[]))] <- 1

my.isoMDS <- isoMDS(dist.dist,k=3)
my.axes <- data.frame(my.isoMDS$points)

#rescale to (roughly) 0 to 1
my.axes[which(is.na(my.axes))] <- -1

rescale <- function(in_vec, new_min = 0, new_max = 1) {
  old_min <- min(in_vec)
  old_max <- max(in_vec)
  out_val <- in_vec - (old_min - new_min)
  out_val <- (out_val / ((old_max - old_min) * new_max))
  return(out_val)
}

my.axes_r <- my.axes
for (col in 1:3) {
  my.axes_r[,col] <- rescale(my.axes_r[,col])
}

my.cols = rgb(my.axes_r)

windows(5,10)
main.text <- paste("Lineage turnover \nresolution:",resolution )
plot(models_grouped$long_group,models_grouped$lat_group,col=my.cols,pch=15,cex=0.8,xlab="longitude",ylab="latitude", main=main.text)

# add a coastline to the map
library(raster)

border.shp <- shapefile(border)
plot(border.shp,add=T)