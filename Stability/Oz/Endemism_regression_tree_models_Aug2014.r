
rm(list=ls())

library(SDMTools)
library(raster)
library(dismo)

colour_range <- function(values, class_count = 20, digits = 3) {
  library(classInt)
  source("~/FromYale/FromTuraco/PhyloSpatial/MammalData/MapFunctions.r")
  
  # define class breaks as quantile
  my.class.fr <- classIntervals(values,n=class_count,style="quantile", digits)
  my.class.fr[[2]] <- round(my.class.fr[[2]],digits)
  
  #set initial colours
  my.pal<-c("light grey","yellow","red") # choose colors
  my.col.fr<-findColours(my.class.fr,my.pal,between="to") # ramp colors based on classInts
  
  return(my.col.fr)
}

base_path     <- 'C:/Users/u3579238/Work/Refugia/'
results.dir   <- paste(base_path,'Results/',sep='')

regions <- data.frame()
do_logistic <- FALSE
logistic_threshold <- 0.95  # for a logistic model to predict membership of top diversity

# whether to do delta AIC table and glmulti
do_glmulti <- FALSE
do_plots   <- FALSE

# whether to fit boosted regression tree models
do_BRT <- TRUE
do_BRT_quantile_plots <- TRUE

i <- 1
regions[i,"region"]        <- 'ALL'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/ALL_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.topo.buf200km/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.topo.buf200km/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output.topo.buf200km/000.asc',sep='')

i <- 2
regions[i,"region"]        <- 'CYP'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.topo.buf200km/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.topo.buf200km/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output.topo.buf200km/000.asc',sep='')

i <- 3
regions[i,"region"]        <- 'AWT'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.topo.buf200km/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.topo.buf200km/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output.topo.buf200km/000.asc',sep='')

i <- 4
regions[i,"region"]        <- 'MEQ'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.topo.buf200km/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.topo.buf200km/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output.topo.buf200km/000.asc',sep='')

i <- 5
regions[i,"region"]        <- 'CEC'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.topo.buf200km/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.topo.buf200km/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output.topo.buf200km/000.asc',sep='')


## EXCLUDING THE SEA REGUION DUE TO INSUFFICIENT DATA
i <- 6
regions[i,"region"]        <- 'SEA'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.topo.buf200km/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.topo.buf200km/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output.topo.buf200km/000.asc',sep='')

i <- 7
regions[i,"region"]        <- 'ALL_not_SEA'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/ALL_rainforest_NOT_SEA.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[1,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.topo.buf200km/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.topo.buf200km/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output.topo.buf200km/000.asc',sep='')

diversities <- data.frame(taxon="lizardfrog",
                          level="lineage", 
                          metric="endemism", 
                          grid="reptfrog_end_lin_18aug14_thresh_01.asc",
                          transform="log",
                          stringsAsFactors = F)
j <- 2
diversities[j,1:4] <- c("lizardfrog","lineage","richness","reptfrog_rich_lin_18aug14_thresh_01.asc")
j <- 3
diversities[j,1:4] <- c("reptile","species","endemism","rept_end_sp.asc")
j <- 4
diversities[j,1:4] <- c("reptile","species","richness","rept_rich_sp.asc")
j <- 5
diversities[j,1:4] <- c("frog","species","endemism","frog_end_sp.asc")
j <- 6
diversities[j,1:4] <- c("frog","species","richness","frog_rich_sp.asc")
#note, lineage richness = species richness.  But it gives richness for only the species with lineages

#output <- data.frame()
k  <- 0

# create a data frame to hold the model results
result_frame <- data.frame(region="",response_description="",predictor_description="",n=0,glm_r2=0,glm_aic=0,glm_delta_aic=0,sarlm_aic=0,sarlm_delta_aic=0,glm_logistic_r2=0,glm_logistic_aic=0,glm_logistic_delta_aic=0,sarlm_logistic_aic=0,sarlm_logistic_delta_aic=0,stringsAsFactors = F)

# SET UP A WINDOW FOR PLOTS
if (do_BRT_quantile_plots) {
  windows(20,12)
  par(mfrow=c(2,4))
}

#Loop through each region, and within it, each diversity metric
#for (i in 1:nrow(regions)) {
#for (i in c(1,7,2,3,4,5,6)) {
for (i in 1) {
  
  # for the all region model, assign region sizes to cells

  #define some basic data
  rf.ras <- raster(regions$veg_grid[i])# read in the vegetation grid
  rf.ras[which(is.finite(rf.ras[]) & rf.ras[] !=1)] <- 0  #set all veg != 1 (rainforests) to 0
  #rf_region.asc     <- read.asc(regions$region[i])
  stabil_static.ras <- raster(regions$stabil_static[i])
  stabil_10m.asc    <- read.asc(regions$stabil_10m[i])
  now_mod.asc       <- read.asc(regions$now_mod[i])
  
  #crop the rainforest raster to match the stability raster
  rf.ras = crop(x=rf.ras,y=stabil_static.ras)
  rf.asc = asc.from.raster(rf.ras)
  stabil_static.asc = asc.from.raster(stabil_static.ras)
  
  j=1 # for now only do lineage endemism
  #for (j in 1:nrow(diversities)) {
  
    #load the diversity result at 0.01 degree resolution and resample to match stability
    div.ras     <-  raster(paste(results.dir,diversities$grid[j],sep=''))
    div_resample.ras <- resample(div.ras,stabil_static.ras,method="bilinear")
    div_resample.asc <- asc.from.raster(div_resample.ras)
    
    #create pos_rf the table of grid cell values for rainforest
    pos_rf <- as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data

    xy_rf <- getXYcoords(rf.asc)
    pos_rf$long  <- xy_rf$x[pos_rf$row]
    pos_rf$lat  <- xy_rf$y[pos_rf$col]
    rm(xy_rf)
  
    pos_rf$rf = rf.asc[cbind(pos_rf$row,pos_rf$col)]            #append the vegetation data
    pos_rf = pos_rf[which(pos_rf$rf==1),]  # filter to rainforest areas
    pos_rf$div = div_resample.asc[cbind(pos_rf$row,pos_rf$col)] #append the diversity data
    pos_rf = pos_rf[which(pos_rf$div > 0),]  # filter to areas with a diversity score
    
    # add the predictors based on rainforest niche models
    pos_rf$stabil_static = stabil_static.asc[cbind(pos_rf$row,pos_rf$col)] #append the static stability scores
    pos_rf$stabil_10m = stabil_10m.asc[cbind(pos_rf$row,pos_rf$col)] #append the dynamic stability scores
    pos_rf$now_mod = now_mod.asc[cbind(pos_rf$row,pos_rf$col)] #append the current rainforest SDM
  
    # add other predictors
    predictors <- data.frame(name="bio1",description="mean annual temperature",path="C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_01.asc",resample=FALSE,stringsAsFactors = F)
    predictors[2,] <- c("bio12","annual precipitation","C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_12.asc",FALSE)
    predictors[3,] <- c("bio17","dry quarter precipitation","C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_17.asc",FALSE)
    predictors[4,] <- c("roughness","topographic complexity (SD of elevation from 9 second DEM)","C:/Users/u3579238/Work/Refugia/Results/elev_sd.asc",FALSE)
    predictors[5,] <- c("region_area","number of pixels in each rainforest region","C:/Users/u3579238/Work/Refugia/Results/rainforest_region_ranges.asc",FALSE)
    # patch size
    
    for (p in 1:nrow(predictors) ) {
      path <- predictors$path[p]
      env.ras <- raster(path)
      if (predictors$resample[p]) {
        env.ras <- resample(env.ras,stabil_static.ras,method="bilinear")
        cat(predictors[p,"name"], "done\n")
      }
      
      pos_rf[,predictors$name[p]] <- extract(env.ras,pos_rf[,c("long","lat")])
    }

    # transform the response variable if needed (as set in the diversities data frame)
    if (diversities[j,"transform"] == "log") {
      pos_rf$div <- log(pos_rf$div)
    }
  
    # rescale the predictors - KEEP THE COLUMN NUMBERS CORRECT
    #pos_rf[,7:14] <- scale(pos_rf[,7:14])
    
    # add a binary variable for membership of top class
    thresh <- quantile(pos_rf$div,logistic_threshold)
    pos_rf$div_binary <- rep(0,nrow(pos_rf))
    pos_rf$div_binary[pos_rf$div >= thresh] <- 1

    n <- nrow(pos_rf)
  
#   # plot current RF suitability v static stability
#   #windows()
#   colours <- colour_range(pos_rf$div, class_count = 20, digits = 1)
#   plot(pos_rf$now_mod,pos_rf$stabil_10m,cex=0.75, cex.lab=1.3, xlab="Current rainforest SDM", ylab="Paleo-stability (10m/yr)", 
#        cex.axis=1.2, pch=20, col=colours, main = regions$region[i])
#   abline(a=0,b=1,col="dark green",lwd=2)
#   
#   #Prepare legend
#   legtext <- names(attr(colours,"table"))  # declare labels
#   legtext <- substr(legtext,2, nchar(legtext)-1) # delete silly brackets
#   legtext <- sub(","," to ",legtext)
#   legshow <- c(20,18,15,12,9,6,3,1) # show selected classes from 50
#   #legend_title <- "log PE"
#   legcols <- attr(colours,"palette")
#   leg_top <- max(pos_rf$stabil_10m,na.rm=T) * 1.05
#   leg_left <- min(pos_rf$now_mod,na.rm=T) - 0.04
#   legend(leg_left,leg_top,legend=legtext[legshow],fill=legcols[legshow], cex=1,bty="n", title = "", pt.cex = 1.2) # selected legend items

  if (do_BRT) {
    library(gbm)
    
    if (regions$region[i] == "ALL" | regions$region[i] == "ALL_not_SEA") {
      predict_columns <- 7:14
    } else {
      predict_columns <- 7:13
    }

    my.brt <- gbm.step(data=pos_rf, gbm.x = predict_columns, gbm.y = 6,
                      family = "gaussian", tree.complexity = 5,
                      learning.rate = 0.005, bag.fraction = 0.5)
    #my_brt.simple <- gbm.simplify(my.brt,n.drops = 2)

    if (do_BRT_quantile_plots)     {
      the.data <- cbind(pos_rf[,c(6,predict_columns)])
      quantile_points <- seq(0.02,0.98,0.02)
      nq <- length(quantile_points)
      
      # create a data frame for the predictor importance results
      quant_results <- data.frame(cbind(quantile_points,pos_rf[1:nq,predict_columns]), row.names = NULL)
      quant_results[,2:ncol(quant_results)] <- NA
  
      row <- 1
      for (q in quantile_points) {
        my_quantile.brt <- gbm(formula = formula(the.data), data=the.data, distribution = list(name="quantile", alpha = q),)
        cat("\nQuantile:",q,"\n")
        qsum <- summary(my_quantile.brt,plotit=F)
        for (r in 1:nrow(qsum)) {
          colname <- as.character(qsum[r,1])
          quant_results[row,colname] <- qsum[r,2]
        }
        row <- row + 1
      }
  
      # plot predictor importance across the quantile range
      plot(quant_results$quantile_points,quant_results$now_mod,type = "l", xlab="Quantile value", ylab="Importance", main = regions$region[i], ylim=c(0,100), lwd=2)
      lines(quant_results$quantile_points,(quant_results$stabil_static + quant_results$stabil_10m),col="blue", lwd=2)
    }
    
    windows()
    summary(my.brt)
    
#     windows()
#     gbm.plot(my.brt, write.title = F, show.contrib = T)
#     
#     windows()
#     gbm.plot.fits(my.brt, use.factor = T)
#     
#     find.int <- gbm.interactions(my.brt)
#     find.int$interactions
    
    windows()
    title <- paste("Influence of current and past rainforest suitability\n on lineage endemism.   Region:",regions$region[i])
    gbm.perspec(my.brt, 2, 3, z.range = c(-10.5,-4.9), theta = -25, perspective = T, x.label="Dynamic stability 10m/yr", y.label = "Current RF model", z.label = "Lineage endemism", main = title)

    # redo this graph formatted for ms
    title <- paste("Influence of current and past rainforest suitability\n on lineage endemism.")
    pdf("gbm perspective stability current rf endemism.pdf",title=title)
    gbm.perspec(my.brt, 2, 3, z.range = c(-10.5,-4.9), theta = -30, phi=20, x.label="Dynamic stability 10m/yr", y.label = "Current RF model", z.label = "Lineage endemism", main = title, cex.lab=1.2)
    dev.off()

    cat("A place to pause")

  }
  
}
