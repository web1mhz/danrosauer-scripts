# a script to create a bias grid for use with Maxent
# Dan Rosauer - September 2013

library(SDMTools)
library(raster, stringr)
rm(list=ls())

########## Parameters ##########
#inputs
#taxon_name        <- "Carlia"
higher_taxon = "skinks"
base.dir = paste("C:/Users/u3579238/work/Phylofest/Models/",higher_taxon,"/",sep="")
samples.filename = paste("species_sites/AllSkinks_ALA.csv",sep="")

#set up sample locations
samples.dir = c(paste(base.dir,"sequence_sites/",sep=""),paste(base.dir,"species_sites/",sep=""))
samples.filter = c("loc.csv","ALA.csv")
samples.lat.col = c(6,2)
samples.long.col= c(7,3)
sample_locations <- data.frame(cbind(samples.dir,samples.filter))
names(sample_locations) <- c("dir","filter")

resolution=0.01

env.filename <- "c:\\Users\\u3579238\\GISData\\EnvironmentGrids\\AusGDMGrids\\ForMaxent\\bio1.asc" # where the original environment grids are stored
output.filename   <- paste(base.dir,"species_models/bias_files/",higher_taxon,"_bias_grid.asc",sep="")
################################

work.dir       <- paste(base.dir,"species_models/bias_files",sep="")
setwd(work.dir)

template.ras <- raster(env.filename)
bias.ras <- raster(extent(template.ras), nrows=nrow(template.ras), ncols=ncol(template.ras))

#read in additional sample sites
for (i in 1:nrow(sample_locations)) {
  dir = as.character(sample_locations$dir[i])
  filter = as.character(sample_locations$filter[i])
  lat.col = samples.lat.col[i]
  long.col = samples.long.col[i]
  files <- list.files(dir,pattern=filter,full.name=TRUE)    

  first=TRUE
  for (file in files) {
    extra_sites <- read.csv(file)
    extra_sites <- extra_sites[,c(long.col,lat.col)]
    names(extra_sites) <- c("Long","Lat")
    if (first) {
      sample_sites <- extra_sites
    } else {
      sample_sites <- rbind(sample_sites,extra_sites)
    }
    sample_sites <- unique(sample_sites)
    cat("Loaded:",file,"\nSites loaded:",nrow(sample_sites),"\n")
    first=FALSE
  }
  rm(extra_sites)
}

#create a shifted set of sample sites so that where a site lies exactly on the cell border, both adjoining cells are included
#rasterize allocates points on the border to the left or lower (or lower-left) cell, hence shifting up and right
shift <- resolution/50 # an arbitrary small shift to put points on the margin into the next cell
shifted_x <- sample_sites[1] + shift
shifted_y <- sample_sites[2] + shift
sample_sites <- rbind(sample_sites,cbind(sample_sites[1],shifted_y),cbind(shifted_x,sample_sites[2]),cbind(shifted_x,shifted_y))
sample_sites <- unique(sample_sites)

bias.ras <- rasterize(x=sample_sites,y=bias.ras,field=1)
writeRaster(bias.ras,output.filename, datatype="INT1S",NAflag=-1,overwrite=TRUE)

cat("\nFinished.\n",date(),"\n")
