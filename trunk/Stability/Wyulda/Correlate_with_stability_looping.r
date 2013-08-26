
rm(list=ls())

library(SDMTools)

base_path     <- 'C:/Users/u3579238/Work/Refugia/'
results.dir   <- paste(base_path,'Results/',sep='')

regions <- data.frame()

i <- 1
regions[i,"region"]        <- 'ALL'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/ALL_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output/000.asc',sep='')

i <- 2
regions[i,"region"]        <- 'CYP'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output/000.asc',sep='')

i <- 3
regions[i,"region"]        <- 'AWT'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.8/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.8/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output8/000.asc',sep='')

i <- 4
regions[i,"region"]        <- 'MEQ'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output/000.asc',sep='')

i <- 5
regions[i,"region"]        <- 'CEC'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability1/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability1/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output1/000.asc',sep='')

i <- 6
regions[i,"region"]        <- 'SEA'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability/shift.10.asc',sep='')
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
j <- 5
diversities[j,1:4] <- c("frog","species","endemism","frog_end_sp.asc")
j <- 6
diversities[j,1:4] <- c("frog","species","richness","frog_rich_sp.asc")
#note, lineage richness = species richness.  But it gives richness for only the species with lineages

output <- data.frame()
k      <- 0

#Loop through each region, and within it, each diversity metric
for (i in 1:nrow(regions)) {

  #define some basic data
  rf.asc = read.asc(regions$veg_grid[i])                 # read in the vegetation grid
  rf.asc[which(is.finite(rf.asc) & rf.asc!=1)] <- 0  #set all veg != 1 (rainforests) to 0
  #rf_region.asc     <- read.asc(regions$region[i])
  stabil_static.asc <- read.asc(regions$stabil_static[i])
  stabil_10m.asc    <- read.asc(regions$stabil_10m[i])
  now_mod.asc       <- read.asc(regions$now_mod[i])
  
  #crop the rainforest raster to match the stability raster
  rf.ras = raster.from.asc(rf.asc)
  stabil_static.ras = raster.from.asc(stabil_static.asc)
  stabil.extent=extent(stabil_static.ras)
  rf.crop.ras = crop(x=rf.ras,stabil.extent)
  rf.asc = asc.from.raster(rf.crop.ras)
  stabil_static.asc = asc.from.raster(stabil_static.ras)
  
  for (j in 1:nrow(diversities)) {
    
    #load the diversity result at 0.01 degree resolution and resample to match stability
    div.asc     <-  read.asc(paste(results.dir,diversities$grid[j],sep=''))
    div.ras     <-  raster.from.asc(div.asc)
    div_resample.ras <- resample(div.ras,stabil_static.ras,method="bilinear")
    div_resample.asc <- asc.from.raster(div_resample.ras)
    
    pos_rf = as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data
    pos_rf$rf = rf.asc[cbind(pos_rf$row,pos_rf$col)]            #append the vegetation data
    pos_rf$stabil_static = stabil_static.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$stabil_10m = stabil_10m.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$now_mod = now_mod.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$div = div_resample.asc[cbind(pos_rf$row,pos_rf$col)] #append the diversity data
    pos_rf = pos_rf[which(pos_rf$rf==1),]  # filter to rainforest areas
    pos_rf = pos_rf[which(pos_rf$div > 0),]  # filter to areas with an endemism score
    
    lm_stabil_static_div  <- lm(pos_rf$div ~ pos_rf$stabil_static)
    cat("\nRegion:", regions$region[i],"\nStability metric: static\n")
    print(diversities[j,1:3])
    print(summary(lm_stabil_static_div))
    # get p value for model from F statistic
    f <- summary(lm_stabil_static_div)$fstatistic
    p_static <- pf(f[1], f[2], f[3], lower=FALSE)
        
    lm_stabil_10m_div     <- lm(pos_rf$div ~ pos_rf$stabil_10m)
    cat("\nRegion:",regions[i,1],"\nStability metric: 10m/yr\n")
    print(diversities[j,1:3])
    print(summary(lm_stabil_10m_div))
    # get p value for model from F statistic
    f <- summary(lm_stabil_10m_div)$fstatistic
    p_10m <- pf(f[1], f[2], f[3], lower=FALSE)    
    
    # see what stability adds above current suitability
    lm_current <- lm(pos_rf$div ~ pos_rf$now_mod)
    current_r2 <- summary(lm_current)$adj.r.squared
    
    lm_cur_stabil_static <- lm(pos_rf$div ~ pos_rf$now_mod + pos_rf$stabil_static)
    cur_stabil_static_r2 <- summary(lm_cur_stabil_static)$adj.r.squared
    extra_r2_static      <- cur_stabil_static_r2 - current_r2
    
    lm_cur_stabil_10m <- lm(pos_rf$div ~ pos_rf$now_mod + pos_rf$stabil_10m)
    cur_stabil_10m_r2 <- summary(lm_cur_stabil_10m)$adj.r.squared
    extra_r2_10m      <- cur_stabil_10m_r2 - current_r2    
    
    k <- k+1   
    output[k,"region"] <- regions$region[i]
    for (l in 1:3) {
      output[k,names(diversities[l])] <- diversities[j,l]
    }
    output[k,"static r2"] <- summary(lm_stabil_static_div)$adj.r.squared
    output[k,"static p"] <- p_static
    output[k,"10myr r2"] <- summary(lm_stabil_10m_div)$adj.r.squared
    output[k,"10myr p"] <- p_10m
    output[k,"cur plus static r2"] <- cur_stabil_static_r2
    output[k,"marginal static r2"] <- extra_r2_static
    output[k,"cur plus 10m r2"] <- cur_stabil_10m_r2 
    output[k,"marginal 10m r2"] <- extra_r2_10m

  }
}

write.csv(output,"Stability_end_rich_cor_with_frogs.csv")

# # now plot current suitability v stability (past suitability), coloured by endemism
# library(maptools)
# library(classInt)
# windows()
# class_count <- 12
# my.class <- classIntervals(pos_rf$lin_end,n=class_count,style="quantile", digits=2)
# my.class_breaks <- round(my.class[[2]],4)
# my.pal <- c("darkblue","green2","yellow","red")
# my.col <-findColours(my.class,my.pal)
# legend_cols <- attr(my.col,"palette")
# plot(pos_rf$now_mod,pos_rf$stabil_static,xlab="Current suitability",ylab="Stability",col=my.col)
# legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)
# 
# # now plot static stability v shifting stability 10m, coloured by endemism
# library(maptools)
# library(classInt)
# windows()
# class_count <- 12
# my.class <- classIntervals(pos_rf$lin_end,n=class_count,style="quantile", digits=2)
# my.class_breaks <- round(my.class[[2]],4)
# my.pal <- c("darkblue","green2","yellow","red")
# my.col <-findColours(my.class,my.pal)
# legend_cols <- attr(my.col,"palette")
# plot(pos_rf$stabil_static,pos_rf$stabil_10m,xlab="Static suitability",ylab="Stability 10m/yr",col=my.col)
# abline(0,1, lwd=2)
# legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)

rm(rf.ras, rf.crop.ras,lin_end_01.asc,lin_end_01.ras)

