

################################################################################
library(SDMTools)

################################################################################
################################################################################

#define directories
# work.dir = 'C:/Users/u3579238/Work/Refugia/Stability/ALL_RF/'; setwd(work.dir)
# veg_grid = 'C:/Users/u3579238/Work/Refugia/Stability/NVIS/ALL_rainforest.asc'
# stabil_static = 'C:/Users/u3579238/Work/Refugia/Stability/ALL_RF/stability/static.sum.cost.asc'
# stabil_10m = 'C:/Users/u3579238/Work/Refugia/Stability/ALL_RF/stability/shift.10.asc'
# now_mod  = 'C:/Users/u3579238/Work/Refugia/Stability/ALL_RF/maxent.output/000.asc'
# lin_end_grid = 'C:/Users/u3579238/Work/Refugia/Results/rept_end_lin_25Jan_thresh_01.asc'

base_path     <- 'C:/Users/u3579238/Work/Refugia/'
results.dir   <- paste(base_path,'Results/',sep='')

regions <- data.frame()

i <- 1
regions[i,"region"]        <- 'MEQ'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/stability/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'_RF/stability/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output/000.asc',sep='')

diversities <- data.frame()
j <- 1
diversities[j,"taxon"]       <- "lizard"
diversities[j,"level"]       <- "lineage"
diversities[j,"metric"]      <- "endemism"
diversities[j,"grid"]        <- "rept_end_lin_25Jan_thresh_01.asc"
j <- 2
diversities[j,1:4] <- c("lizard","lineage","richness","rept_rich_lin_25Jan_thresh_01.asc")
j <- 3
diversities[j,1:4] <- c("reptile","species","endemism","rept_end_sp.asc")
j <- 4
diversities[j,1:4] <- c("reptile","species","richness","rept_rich_sp.asc")
#note, lineage richness = species richness.  But it gives richness for only the species with lineages

#Loop through each region, and within it, each diversity metric
for (i in 1:nrow(regions) {

  #define some basic data
  rf.asc = read.asc(regions$veg_grid[i]])                 # read in the vegetation grid
  rf.asc[which(is.finite(rf.asc) & rf.asc!=1)] = 0  #set all veg != 1 (rainforests) to 0
  rf_region.asc     <- read.asc(regions$region[i])
  stabil_static.asc <- read.asc(regions$stabil_static[i})
  stabil_10m.asc    <- read.asc(regions$stabil_10m[i])
  now_mod.asc       <- read.asc(regions$now_mod[i])
  
  #crop the rainforest raster to match the stability raster
  rf.ras = raster.from.asc(rf.asc)
  stabil_static.ras = raster.from.asc(stabil_static.asc)
  stabil.extent=extent(stabil_static.ras)
  rf.crop.ras = crop(x=rf.ras,stabil.extent)
  rf.asc = asc.from.raster(rf.crop.ras)
  stabil_static.asc = asc.from.raster(stabil_static.ras)
  
  #crop the region raster to match the stability raster
  rf_region.ras = raster.from.asc(rf_region.asc)
  rf_region.crop.ras = crop(x=rf_region.ras,stabil.extent)
  rf_region.asc = asc.from.raster(rf_region.crop.ras)
  
  for (j in 1:nrow(diversities)) {
    
    #load the diversity result at 0.01 degree resolution and resample to match stability
    div.asc     <-  read.asc(diversities$grid)
    div.ras     <-  raster.from.asc(div.asc)
    div_resample.ras <- resample(div.ras,stabil_static.ras,method="bilinear")
    div_resample.asc <- asc.from.raster(div_resample.ras)
    
    pos_rf = as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data
    pos_rf$rf = rf.asc[cbind(pos_rf$row,pos_rf$col)]            #append the vegetation data
    pos_rf$region = rf_region.asc[cbind(pos_rf$row,pos_rf$col)] #append the region data
    pos_rf$stabil_static = stabil_static.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$stabil_10m = stabil_10m.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$now_mod = now_mod.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$div = div_resample.asc[cbind(pos_rf$row,pos_rf$col)] #append the diversity data
    pos_rf = pos_rf[which(pos_rf$rf==1),]  # filter to rainforest areas
    pos_rf = pos_rf[which(pos_rf$div > 0),]  # filter to areas with an endemism score
    
  }

}
# 
# #load the lineage richness result at 0.01 degree resolution and resample to match stability
# lin_rich_01.asc = read.asc(lin_rich_grid)
# lin_rich_01.ras = raster.from.asc(lin_rich_01.asc)
# lin_rich_resample.ras = resample(lin_rich_01.ras,stabil_static.ras,method="bilinear")
# lin_rich_resample.asc = asc.from.raster(lin_rich_resample.ras)

#load the species endemism result at 0.01 degree resolution and resample to match stability
# spec_end.asc = read.asc(sp)
# lin_end_01.ras = raster.from.asc(lin_end_01.asc)
# lin_end_resample.ras = resample(lin_end_01.ras,stabil_static.ras,method="bilinear")
# lin_end_resample.asc = asc.from.raster(lin_end_resample.ras)
# 
lm_bio_stabil <- lm(pos_rf$lin_end ~ pos_rf$stabil_static)
print(summary(lm_bio_stabil))

lm_bio_stabil <- lm(pos_rf$lin_end ~ pos_rf$stabil_10m)
print(summary(lm_bio_stabil))

lm_bio_stabil <- lm(pos_rf$lin_rich ~ pos_rf$stabil_static)
print(summary(lm_bio_stabil))

lm_bio_stabil <- lm(pos_rf$lin_rich ~ pos_rf$stabil_10m)
print(summary(lm_bio_stabil))

# see what stability adds above current suitability
lm_bio_stabil <- lm(pos_rf$lin_end ~ pos_rf$now_mod)
print(summary(lm_bio_stabil))

lm_bio_stabil <- lm(pos_rf$lin_end ~ pos_rf$now_mod + pos_rf$stabil_static)
print(summary(lm_bio_stabil))

lm_bio_stabil <- lm(pos_rf$lin_end ~ pos_rf$now_mod + pos_rf$stabil_10m)
print(summary(lm_bio_stabil))


# now plot current suitability v stability (past suitability), coloured by endemism
library(maptools)
library(classInt)
windows()
class_count <- 12
my.class <- classIntervals(pos_rf$lin_end,n=class_count,style="quantile", digits=2)
my.class_breaks <- round(my.class[[2]],4)
my.pal <- c("darkblue","green2","yellow","red")
my.col <-findColours(my.class,my.pal)
legend_cols <- attr(my.col,"palette")
plot(pos_rf$now_mod,pos_rf$stabil_static,xlab="Current suitability",ylab="Stability",col=my.col)
legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)

# now plot static stability v shifting stability 10m, coloured by endemism
library(maptools)
library(classInt)
windows()
class_count <- 12
my.class <- classIntervals(pos_rf$lin_end,n=class_count,style="quantile", digits=2)
my.class_breaks <- round(my.class[[2]],4)
my.pal <- c("darkblue","green2","yellow","red")
my.col <-findColours(my.class,my.pal)
legend_cols <- attr(my.col,"palette")
plot(pos_rf$stabil_static,pos_rf$stabil_10m,xlab="Static suitability",ylab="Stability 10m/yr",col=my.col)
abline(0,1, lwd=2)
legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)


# relate lineage endemism to species richness (lineage richness = species richness, but used here to limit to species with lineage data)
lm_lin_sp


rm(rf.ras, rf.crop.ras,lin_end_01.asc,lin_end_01.ras)

# NEXT REDO FOR EAST REGION SEPARATELY
# AND WITH SINGLE REGION STABILITY MODEL
# AND WITH OTHER ENDEMISM MEASURES (FROG, REPT_LIN)

