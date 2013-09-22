
rm(list=ls())

library(SDMTools)
library(raster)

# parameters #############################################################
input.dir   = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models_aligned_rf_strict/'; setwd(input.dir)
output.dir  = 'C:/Users/u3579238/Work/Refugia/Results/outputs_11Jun'
template_ext ='C:/Users/u3579238/GISData/EnvironmentGrids/AusGDMGrids/ForMaxent/bio1.asc'
base_path     <- 'C:/Users/u3579238/Work/Refugia/'
results.dir   <- paste(base_path,'Results/',sep='')
file.pattern  <- '.+asc$'  #regex
threshold     <- 0.9

region_lat_boundaries <- c(-9.7,-14.85,-19.6,-23.4,-32.44, -44.3)
##########################################################################

regions <- data.frame()

i <- 1
regions[i,"region"]        <- 'ALL'
regions[i,"veg_grid_coarse"]      <- paste(base_path,'Stability/NVIS/ALL_rainforest.asc',sep='')
regions[i,"veg_grid_fine"] <- paste(base_path,'Stability/NVIS/nvis4_1_aust_mvs_pre_geo01_rf_only.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output/000.asc',sep='')

diversities <- data.frame()
j <- 1
diversities[j,"taxon"]       <- "lizard"
diversities[j,"level"]       <- "lineage"
diversities[j,"metric"]      <- "endemism"
diversities[j,"grid"]        <- "rept_end_lin_11Jun_thresh_01.asc"
j <- 2; diversities[j,1:4] <- c("lizard","lineage","richness","rept_rich_lin_25Jan_thresh_01.asc")
j <- 3; diversities[j,1:4] <- c("reptile","species","endemism","rept_end_sp.asc")
j <- 4; diversities[j,1:4] <- c("reptile","species","richness","rept_rich_sp.asc")
j <- 5; diversities[j,1:4] <- c("frog","species","endemism","frog_end_sp.asc")
j <- 6; diversities[j,1:4] <- c("frog","species","richness","frog_rich_sp.asc") 

# import the lineage models at 0.01 degrees and clip them to RF
i <- 1
# list the lineage models
input_files= list.files(path = input.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE,
                        ignore.case = TRUE, include.dirs = FALSE)
#get the extent of the models and crop the veg grid to match
filepath=paste(input.dir,input_files[1],sep='')
lin_mod.ras <- raster(filepath)
rf.ras      <- raster(regions$veg_grid_fine[i])                 # read in the vegetation grid
rf.ras      <- crop(rf.ras,lin_mod.ras)
rf.asc      <- asc.from.raster(rf.ras)
pos = as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE))         #get all points that have data
pos$rf = rf.asc[cbind(pos$row,pos$col)]                             #append the vegetation data
pos$rf[pos$rf>0] <- 1

#add lat, long
latlong = getXYcoords(rf.asc)
pos$lat <- latlong$y[pos$col]
pos$long <- latlong$x[pos$row]

# add in the coarse grain data
rf_coarse.asc     <- read.asc(regions$veg_grid_coarse[i])                 # read in the vegetation grid
stabil_static.asc <- read.asc(regions$stabil_static[i])
stabil_static.ras <- raster.from.asc(stabil_static.asc)
stabil_10m.asc    <- read.asc(regions$stabil_10m[i])
now_mod.asc       <- read.asc(regions$now_mod[i])

pos$rf_coarse  <- extract.data(pos[c("long","lat")],rf_coarse.asc)
pos$now_mod  <- extract.data(pos[c("long","lat")],now_mod.asc)
pos$stabil_static  <- extract.data(pos[c("long","lat")],stabil_static.asc)
pos$stabil_10m  <- extract.data(pos[c("long","lat")],stabil_10m.asc)

# identify each pixel by region based on latitude
for (region in 1:(length(region_lat_boundaries)-1)) {
  lat_max = region_lat_boundaries[region]
  lat_min = region_lat_boundaries[region+1]
  pos$region[pos$lat <= lat_max & pos$lat >= lat_min] <- region
}

first_lin_col = ncol(pos) + 1

# read in each of the lineage models
for (tfile in input_files) {
  filepath=paste(input.dir,tfile,sep='')
  dataname = unlist(strsplit(tfile,".asc")[1])[1]  #define the column name
  
  #load the lineage model at 0.01 degree resolution and make the values a column. Each cell is a row
  lin_mod.asc <- read.asc(filepath)
  pos[,dataname] <- lin_mod.asc[cbind(pos$row,pos$col)] #append the lineage model data
  cat(dataname,"\tdone\n")
}

# grid cell stats
cols = ncol(pos)
pos$total = apply(pos[first_lin_col:cols],1,'sum',na.rm=T)      # analogous to richness
lineage_sums = apply(pos[first_lin_col:cols],2,'sum',na.rm=T)

# the endemism calculation, cumulative model values and presence/absence
pos.end = pos[1:cols]
pos.cum = pos.end
pos.occ = pos.end
for (lineage in names(pos[first_lin_col:cols])) {
  if (lineage_sums[lineage] > 0) {    # avoid division by zero if a species column sums to 0
    pos.end[lineage] <- (pos.end[lineage] / lineage_sums[lineage])
    pos.end[is.na(pos.end[lineage]),lineage] <- 0
    order.decr <- order(pos.end[lineage],decreasing=TRUE)
    
    #cumulative model values
    pos.cum[order.decr,lineage] <- cumsum(pos.end[order.decr,lineage])
    
    #presence / absence values
    pos.occ[pos.cum[lineage] <= threshold,lineage] <- 1
    pos.occ[pos.cum[lineage] > threshold,lineage] <- 0
  }
}

pos$we = apply(pos.end[first_lin_col:cols],1,'sum',na.rm=T)
lineage_ranges = apply(pos.occ[first_lin_col:cols],2,'sum',na.rm=T)
lineage_ranges = as.data.frame(lineage_ranges)
lineage_ranges$lineage <- rownames(lineage_ranges)
names(lineage_ranges)[1] <- "range"
lineage_ranges <- data.frame(lineage_ranges[2:1],row.names=NULL)

region_area <- 1:5

#divide lineage ranges by region
for (region in 1:(length(region_lat_boundaries)-1)) {
  lineage_ranges[region+2] <- rep(0,nrow(lineage_ranges))
  names(lineage_ranges)[region+2] <- paste("region_",region,sep="")
  pos.region <- pos.occ[pos.occ$region == region,]
  region_area[region] <- length(pos.region[,1])
  for (lineage in lineage_ranges$lineage) {
    lineage_ranges[lineage_ranges$lineage==lineage,region+2] <- sum(pos.region[,lineage],na.rm=T)
  }
}

#calculate proportion of range in each region - for presence-absence
lineage_ranges_proportion <- lineage_ranges[,3:7]
lineage_ranges_proportion <- round(lineage_ranges_proportion / lineage_ranges$range,2)
lineage_ranges_proportion <- cbind(lineage_ranges[,1:2],lineage_ranges_proportion)
lineage_ranges_proportion$endemic <- apply(lineage_ranges_proportion[,3:7],1,max)==1
#lineage_ranges_proportion[lineage_ranges_proportion$endemic==TRUE,"end_region"] <- 
for (row in 1:nrow(lineage_ranges_proportion)) {
  region <- which(lineage_ranges_proportion[row,3:7]==1)
  if (length(region) > 0) {
    lineage_ranges_proportion$end_region[row] <- region
  } else {
    lineage_ranges_proportion$end_region[row] <- 0    
  }
}

#divide lineage **likelihoods** by region
for (region in 1:(length(region_lat_boundaries)-1)) {
  lineage_ranges[region+2] <- rep(0,nrow(lineage_ranges))
  names(lineage_ranges)[region+2] <- paste("region_",region,sep="")
  pos.region <- pos.occ[pos.occ$region == region,]
  region_area[region] <- length(pos.region[,1])
  for (lineage in lineage_ranges$lineage) {
    lineage_ranges[lineage_ranges$lineage==lineage,region+2] <- sum(pos.region[,lineage],na.rm=T)
  }
}

#calculate proportion of range in each region - for **model likelihood**
lineage_ranges_proportion <- lineage_ranges[,3:7]
lineage_ranges_proportion <- round(lineage_ranges_proportion / lineage_ranges$range,2)
lineage_ranges_proportion <- cbind(lineage_ranges[,1:2],lineage_ranges_proportion)
lineage_ranges_proportion$endemic <- apply(lineage_ranges_proportion[,3:7],1,max)==1
lineage_ranges_proportion[lineage_ranges_proportion$endemic==TRUE,"end_region"] <- 
for (row in 1:nrow(lineage_ranges_proportion)) {
  region <- which(lineage_ranges_proportion[row,3:7]==1)
  if (length(region) > 0) {
    lineage_ranges_proportion$end_region[row] <- region
  } else {
    lineage_ranges_proportion$end_region[row] <- 0    
  }
}

names(lineage_ranges_proportion)[3:7] <- paste(names(lineage_ranges_proportion)[3:7],"_prop",sep="")
lineage_ranges <- cbind(lineage_ranges,lineage_ranges_proportion[,3:ncol(lineage_ranges_proportion)])

pos.cum[first_lin_col:cols] <- round(pos.cum[first_lin_col:cols],5)

# output fine-grain results
orig.dir=getwd()
setwd(output.dir)
write.csv(pos,"rf_lin_models.csv")
write.csv(pos.cum,"rf_lin_models_cumulative.csv")
write.csv(pos.occ,paste("rf_lin_models_occ_",threshold,".csv",sep=""))
write.csv(as.data.frame(lineage_ranges),paste("lineage_ranges_",threshold,".csv",sep=""))
write.csv(as.data.frame(region_area),"region_areas.csv")
setwd(orig.dir)
