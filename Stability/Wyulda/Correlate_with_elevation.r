
rm(list=ls())

library(SDMTools)

base_path     <- 'C:/Users/u3579238/Work/Refugia/'
results.dir   <- paste(base_path,'Results/',sep='')

regions <- data.frame()

i <- 1
regions[i,"region"]        <- 'CEC'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"dem"]           <- 'C:/Users/u3579238/GISData/EnvironmentGrids/Topo/gtopo_1km.asc'
#regions[i,"dem"]           <- 'C:/Users/u3579238/Work/Refugia/Results/below_max_elev_015.asc'
#regions[i,"dem"]           <- 'C:/Users/u3579238/Work/Refugia/Results/bel_loc_max_015_ab_min_250.asc'

diversities <- data.frame()
j <- 1
diversities[j,"taxon"]       <- "lizard"
diversities[j,"level"]       <- "lineage"
diversities[j,"metric"]      <- "endemism"
diversities[j,"grid"]        <- "rept_end_lin_25Jan_thresh_01.asc"
j <- 2
diversities[j,1:4] <- c("frog","species","endemism","frog_end_sp.asc")
# j <- 3
# diversities[j,1:4] <- c("reptile","species","endemism","rept_end_sp.asc")
# j <- 4
# diversities[j,1:4] <- c("reptile","species","richness","rept_rich_sp.asc")
#note, lineage richness = species richness.  But it gives richness for only the species with lineages

output <- data.frame()
k      <- 0

#Loop through each region, and within it, each diversity metric
#for (i in 1:nrow(regions)) {

  #define some basic data
  rf.asc = read.asc(regions$veg_grid[i])                 # read in the vegetation grid
  rf.asc[which(is.finite(rf.asc) & rf.asc!=1)] <- 0  #set all veg != 1 (rainforests) to 0
  dem.asc          <- read.asc(regions$dem[i])  
  
  #define the analysis extent
  dem.ras = raster.from.asc(dem.asc)
  dem.extent=extent(dem.ras)
  rf.ras = raster.from.asc(rf.asc) 
  rf.extent = extent(rf.ras)
  new.extent = rf.extent
#   
#   #CEC
#   new.extent@xmin <- 151
#   new.extent@xmax <- 153.8  
#   new.extent@ymin <- -32.6
#   new.extent@ymax <- -23.9  

#   #Border Conondale CEC
#   new.extent@xmin <- 151.4
#   new.extent@xmax <- 153.8  
#   new.extent@ymin <- -28.8
#   new.extent@ymax <- -26.4  

  #Border CEC
  new.extent@xmin <- 152.0
  new.extent@xmax <- 153.8  
  new.extent@ymin <- -28.8
  new.extent@ymax <- -27.7  
# 
#   #AWT
#   new.extent@xmin <- 144.8
#   new.extent@xmax <- 147.4  
#   new.extent@ymin <- -19.7
#   new.extent@ymax <- -14.7 
# 
#   #MEQ
#   new.extent@xmin <- 147.6
#   new.extent@xmax <- 149.6  
#   new.extent@ymin <- -22.2
#   new.extent@ymax <- --19.9    

  #crop the rainforest and dem rasters to match new extent
  rf.crop.ras = crop(x=rf.ras,new.extent)
  dem.crop.ras = crop(x=dem.ras,new.extent)
  
  # resample the dem to match the diversity layers
  div.asc     <-  read.asc(paste(results.dir,diversities$grid[1],sep=''))  
  div.ras     <-  raster.from.asc(div.asc) 
  div.crop.ras <- crop(x=div.ras,new.extent)
  dem.resample.ras <- resample(dem.crop.ras,div.crop.ras,method="bilinear")
  dem.asc <- asc.from.raster(dem.resample.ras)

  # resample the rf to match the diversity layers
  rf.resample.ras <- resample(rf.crop.ras,div.crop.ras,method="ngb")
  rf.asc   <- asc.from.raster(rf.resample.ras)
  
  #j=1  #temp
  for (j in 1:nrow(diversities)) {
    div.asc     <-  read.asc(paste(results.dir,diversities$grid[j],sep=''))  
    div.ras     <-  raster.from.asc(div.asc) 
    div.crop.ras <- crop(x=div.ras,new.extent)
    div.asc <- asc.from.raster(div.crop.ras)
    
    pos_dem = as.data.frame(which(is.finite(dem.asc),arr.ind=TRUE)) #get all points that have data
    pos_dem$rf = rf.asc[cbind(pos_dem$row,pos_dem$col)]            #append the vegetation data
    pos_dem$dem = dem.asc[cbind(pos_dem$row,pos_dem$col)] #append the dem data
    pos_dem$div = div.asc[cbind(pos_dem$row,pos_dem$col)] #append the diversity data
    pos_dem = pos_dem[which(pos_dem$rf==1),]  # filter to rainforest areas
    #pos_dem = pos_dem[which(pos_dem$div > 0.0001),]  # filter to areas with a non-negligible endemism score
    
    #add percentiles
    percentiles <- ecdf(pos_dem$div)
    pos_dem$div_perc <- percentiles(pos_dem$div)
    
    graph_title <- paste("Lineage Endemism v Elevation for\n",regions[i,"region"])
        
#     # plot elevation v div measure
#     windows()
#     plot(pos_dem$dem,pos_dem$div,ylab=paste(diversities[j,1],diversities[j,2],diversities[j,3],sep=" "),xlab="elevation", main = graph_title)
#     
#     # plot elevation v div measure for top 10% of div
#     windows()
#     plot(pos_dem$dem[pos_dem$div_perc>=0.9],pos_dem$div[pos_dem$div_perc>=0.9],ylab=paste(diversities[j,1],diversities[j,2],diversities[j,3],sep=" "),xlab="elevation", main = graph_title)

#     # histogram
#     graph_title <- paste("Elevation of all rainforest and top 25%, 10%\n",regions[i,"region"],diversities[j,1],diversities[j,2],diversities[j,3])
#     windows(6,6)
#     classcount=15
#     firsthist = hist(pos_dem$dem,breaks=classcount,freq=T,xlab="Elevation", col="grey", main=graph_title)
#     hist(pos_dem[pos_dem$div_perc>=0.75,"dem"],freq=T,add=T,col="blue",breaks=firsthist$breaks)    
#     hist(pos_dem[pos_dem$div_perc>=0.90,"dem"],freq=T,add=T,col="green",breaks=firsthist$breaks)    
#     
    # box and whiskers
    graph_title <- paste("Elevation of rainforest by score\n",regions[i,"region"],diversities[j,1],diversities[j,2],diversities[j,3])
    windows(6,6)
    elev_end_75 <- pos_dem$dem[(pos_dem$div_perc>=0.75) & (pos_dem$div_perc < 0.85)]
    elev_end_85 <- pos_dem$dem[(pos_dem$div_perc>=0.85) & (pos_dem$div_perc < 0.95)]
    elev_end_95 <- pos_dem$dem[(pos_dem$div_perc>=0.95)]    
    boxplot(pos_dem$dem,elev_end_75,elev_end_85,elev_end_95,names=c("All RF","75-85%","85-95%",">95%"),main=graph_title,ylab="Elevation")
    
    
#    abline(v=median(pos_dem$dem,na.rm=T),col="grey",lwd=3)
#    #abline(v=median(pos_dem$dem[pos_dem$div_perc>=0.75],na.rm=T),col="blue",lwd=3)
#    abline(v=median(pos_dem$dem[pos_dem$div_perc>=0.9],na.rm=T),col="green",lwd=3)
#     
#     lm_dem_div  <- lm(pos_dem$div ~ pos_dem$dem)
#     cat("\nRegion:", regions$region[i],"\nMetric: elevation\n")
#     print(diversities[j,1:3])
#     print(summary(lm_dem_div))
#     
#     lm_stabil_10m_div     <- lm(pos_dem$div ~ pos_dem$stabil_10m)
#     cat("\nRegion:",regions[i,1],"\nStability metric: 10m/yr\n")
#     print(diversities[j,1:3])
#     print(summary(lm_stabil_10m_div))
#     
#     # see what stability adds above current suitability
#     lm_current <- lm(pos_dem$div ~ pos_dem$now_mod)
#     current_r2 <- summary(lm_current)$adj.r.squared
#     
#     lm_cur_stabil_static <- lm(pos_dem$div ~ pos_dem$now_mod + pos_dem$stabil_static)
#     cur_stabil_static_r2 <- summary(lm_cur_stabil_static)$adj.r.squared
#     extra_r2_static      <- cur_stabil_static_r2 - current_r2
#     
#     lm_cur_stabil_10m <- lm(pos_dem$div ~ pos_dem$now_mod + pos_dem$stabil_10m)
#     cur_stabil_10m_r2 <- summary(lm_cur_stabil_10m)$adj.r.squared
#     extra_r2_10m      <- cur_stabil_10m_r2 - current_r2    
#     
#     k <- k+1   
#     output[k,"region"] <- regions$region[i]
#     for (l in 1:3) {
#       output[k,names(diversities[l])] <- diversities[j,l]
#     }
#     output[k,"static r2"] <- summary(lm_stabil_static_div)$adj.r.squared
#     output[k,"10myr r2"] <- summary(lm_stabil_10m_div)$adj.r.squared
#     output[k,"cur plus static r2"] <- cur_stabil_static_r2
#     output[k,"marginal static r2"] <- extra_r2_static
#     output[k,"cur plus 10m r2"] <- cur_stabil_10m_r2 
#     output[k,"marginal 10m r2"] <- extra_r2_10m
  }
#}

# # now plot current suitability v stability (past suitability), coloured by endemism
# library(maptools)
# library(classInt)
# windows()
# class_count <- 12
# my.class <- classIntervals(pos_dem$lin_end,n=class_count,style="quantile", digits=2)
# my.class_breaks <- round(my.class[[2]],4)
# my.pal <- c("darkblue","green2","yellow","red")
# my.col <-findColours(my.class,my.pal)
# legend_cols <- attr(my.col,"palette")
# plot(pos_dem$now_mod,pos_dem$stabil_static,xlab="Current suitability",ylab="Stability",col=my.col)
# legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)
# 
# # now plot static stability v shifting stability 10m, coloured by endemism
# library(maptools)
# library(classInt)
# windows()
# class_count <- 12
# my.class <- classIntervals(pos_dem$lin_end,n=class_count,style="quantile", digits=2)
# my.class_breaks <- round(my.class[[2]],4)
# my.pal <- c("darkblue","green2","yellow","red")
# my.col <-findColours(my.class,my.pal)
# legend_cols <- attr(my.col,"palette")
# plot(pos_dem$stabil_static,pos_dem$stabil_10m,xlab="Static suitability",ylab="Stability 10m/yr",col=my.col)
# abline(0,1, lwd=2)
# legend(x="topleft",legend=my.class_breaks[1:class_count+1],fill=legend_cols)

rm(rf.ras, rf.crop.ras,dem.asc,dem.crop.ras,dem.ras,dem.resample.ras,div.ras,div.asc,div.crop.ras,rf.asc,rf.resample.ras)

