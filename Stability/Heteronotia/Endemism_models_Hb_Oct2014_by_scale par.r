
# this script is models environmental correlates of PE in Heteronotia binoei
# October 2014

rm(list=ls())

library(SDMTools)
library(raster)
library(glmulti)
library(dismo)
library(spatial.tools)
library(foreach)

source("C:/Users/u3579238/Work/Software/danrosauer-scripts/Phylofest_code/summarise_to_different_raster.r")

# first some functions - main script follows below
glm_process <- function(result_frame,
                        model_number,
                        predictor_text,
                        response_text,
                        response_text_logistic = NULL, # the response variable name, not values                        the_data,
                        do_spatial_LM = FALSE,
                        weight_list = NULL,
                        do_logistic = FALSE,
                        the_data
)
{
  cat("\n\n** ",model_number,result_frame[model_number,1],result_frame[model_number,2],result_frame[model_number,3]," **\n")

  the.formula <- paste(response_text,"~",predictor_text)
  glm_gauss <- glm(the.formula, data=the_data, family="gaussian")
  print(summary(glm_gauss))
  cat("\nGaussian model\naic:",glm_gauss$aic)
  dev_exp <- (glm_gauss$null.deviance-glm_gauss$deviance)/glm_gauss$null.deviance
  cat("\nDeviance explained:",dev_exp,"\n")
  result_frame[model_number,"glm_r2"] <- round(dev_exp,4)
  result_frame[model_number,"glm_aic"] <- round(glm_gauss$aic,1)

  if (do_logistic) {
    the.formula <- paste(response_text_logistic,"~",predictor_text)
    glm_logistic <- glm(the.formula, data=the_data, family = "binomial")
    print(summary(glm_logistic))
    cat("\nLogistic model \naic:",glm_logistic$aic)
    dev_exp <- (glm_logistic$null.deviance-glm_logistic$deviance)/glm_logistic$null.deviance
    cat("\nDeviance explained:",dev_exp,"\n")
    result_frame[model_number,"glm_logistic_r2"] <- round(dev_exp,4)
    result_frame[model_number,"glm_logistic_aic"] <- round(glm_logistic$aic,1)
  }

  if (do_spatial_LM) {

    SAR_gauss <- errorsarlm(glm_gauss, data=the_data, listw=weight_list, quiet=FALSE, na.omit, zero.policy=TRUE, tol.solve=1e-11)
    #calculate AIC
    SAR_gauss.aic <- (-2*SAR_gauss$LL)+(2*SAR_gauss$parameters)

    print(summary(SAR_gauss))
    cat("\nSAR gaussian \naic:",SAR_gauss.aic)
    result_frame[model_number,"sarlm_aic"] <- round(SAR_gauss.aic,1)

    if (do_logistic) {
      SAR_logistic <- errorsarlm(glm_logistic, data=the_data, listw=weight_list, quiet=FALSE, na.omit, zero.policy=TRUE, tol.solve=1e-11)
      #calculate AIC
      SAR_logistic.aic <- (-2*SAR_logistic$LL)+(2*SAR_logistic$parameters)

      print(summary(SARer_with_exp))
      cat("\nSARLM logistic \naic:",SAR_logistic.aic)
      result_frame[model_number,10] <- round(SAR_logistic.aic,1)
    }

  }

  return(result_frame)
}

# a combine function to use with foreach
my_combine <- function(...) {
  dots <- rbind(...)
  return(dots)
}


base_path     <- 'C:/Users/u3579238/Work/AMT/Models/diversity/H_binoei/'
results.dir   <- paste(base_path,'correlations/',sep='')

resolution    <- 1.75  # resolution for analysis
mask          <- "C:/Users/u3579238/Work/AMT/Maps/mon_trop_poly_clip.shp"

do_logistic <- FALSE
logistic_threshold <- 0.95  # for a logistic model to predict membership of top diversity

#whether to do spatial autocorrelation, and at what radius
do_spatial_LM <- FALSE
spatial_LM_radius <- 100

do_plots   <- FALSE

# whether to fit boosted regression tree models
do_BRT <- FALSE

diversities <- data.frame(taxon="H_binoei",
                          level="lineage",
                          metric="endemism",
                          grid="Hb_PE.asc",
                          transform="log",
                          stringsAsFactors = F)

diversities[2,1:4] <- data.frame(taxon="H_binoei",
                                  level="lineage",
                                  metric="richness",
                                  grid="Hb_SR.asc",
                                  transform="",
                                  stringsAsFactors = F)



# LOOP HERE TO CHECK ACROSS SCALES
scale_results <- data.frame(resolution=0,r2=0,r2_lgm=0)
scale_results <- scale_results[-1,]
i=0

resolutions <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1,1.2,1.4,1.6,1.8,2,2.25,2.5,2.75,3)

#resolutions <- c(1.5,1.75,2)

# create a data frame to hold the model results
glm_blank_results <- data.frame(resolution=0, n=0, r2=0,formula="",stringsAsFactors = F)
# add predictor columns
glm_blank_results <- cbind(glm_blank_results,data.frame(bio1=0,bio12=0,bio17=0,bio1DIFF=0,bio12PROP=0,bio17PROP=0,AET=0,Elev=0,Elev_sd=0, max_fpar_quantile=0, bio12_max_fpar_quantile=0,stringsAsFactors = F))

cpu=5
sfQuickInit(cpus=cpu)

cat("\n\nSTARTING PARALLEL LOOP on",cpu,"processors\n\n")

scale_results <- foreach(i=1:length(resolutions), .combine = my_combine, .packages=c("raster","spatial.tools","SDMTools","glmulti")) %dopar% {
#scale_results <- foreach(i=1:length(resolutions), .packages=c("raster","spatial.tools","SDMTools","glmulti")) %dopar% {

#scale_results <- glm_blank_results[-1,]

#for (i in 1:length(resolutions)) {

  resolution <- resolutions[i]

  glm_one_res_results <- glm_blank_results

  if (do_spatial_LM) {
    library(spdep)
  }

  #define some basic data
  mask.shp <- shapefile(mask)
  the.ext <- extent(mask.shp)
  template.ras   <- raster(ext=the.ext, resolution=resolution)
  template.ras[] <- rep(1,ncell(template.ras))
  template.ras <- mask(template.ras,mask.shp)

  j=1 # for now only do lineage endemism
  #for (j in 1:nrow(diversities)) {

    #load the diversity result at 0.01 degree resolution and resample to match stability
    div.ras     <-  raster(paste(base_path,diversities$grid[j],sep=''))
    div_resample.ras <- resample(div.ras,template.ras,method="bilinear")
    div_resample.ras <- crop(div_resample.ras, template.ras)
    div_resample.asc <- asc.from.raster(div_resample.ras)

    #create pos_df the table of grid cell values for rainforest
    pos_df <- as.data.frame(which(is.finite(div_resample.asc),arr.ind=TRUE)) #get all points that have data

    xy_df <- getXYcoords(div_resample.asc)
    pos_df$long  <- xy_df$x[pos_df$row]
    pos_df$lat  <- xy_df$y[pos_df$col]
    rm(xy_df)

    pos_df$div = div_resample.asc[cbind(pos_df$row,pos_df$col)] #append the diversity data
    pos_df = pos_df[which(pos_df$div > 0),]  # filter to areas with a diversity score

    # add the environment predictors

  predictors <- data.frame(name="bio1",description="mean annual temperature",path="C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_01.asc", resample="standard", stringsAsFactors = F)
  predictors[2,] <- c("bio12","annual precipitation","C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_12.asc", resample="standard")
  predictors[3,] <- c("bio17","dry quarter precipitation","C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_17.asc", resample="standard")
  predictors[4,] <- c("bio1DIFF","change in mean annual temperature since 21ka","C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim_diff/bioclim_01_Now-LGM.asc", resample="standard")
  predictors[5,] <- c("bio12PROP","proportional change in annual precipitation since 21ka","C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim_diff/bioclim_12_log_Now-LGM.asc", resample="standard")
  predictors[6,] <- c("bio17PROP","proportional change in dry quarter precipitation since 21ka","C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim_diff/bioclim_17_log_Now-LGM.asc", resample="standard")
  predictors[7,] <- c("AET","actual evapotranspiration","C:/Users/u3579238/GISData/EnvironmentGrids/AusGDMGrids/flt_grids/rsaet1kcfc.flt", resample="standard")
  predictors[8,] <- c("Elev","elevation","C:/Users/u3579238/GISData/EnvironmentGrids/Topo/gtopo_1km.asc", resample="standard")
  predictors[9,] <- c("Elev","elevation","C:/Users/u3579238/GISData/EnvironmentGrids/Topo/gtopo_1km.asc", resample="sd")
  predictors[10,] <- c("max_fpar","max_fpar","C:/Users/u3579238/GISData/Mackey/Greenspots_Apr13/fgreen00_12_yr_min_max005deg_AMT.asc", resample="quantile")
  # at present, quantile returns the 75th percentile - difficult to change


  # for each raster, identify the values in each template cell, and summarise
  for (p in 1:nrow(predictors) ) {
    path <- predictors$path[p]
    env.ras <- raster(path)
    projection(env.ras) <- projection(template.ras)

    if (predictors$resample[p] == "standard") {
      env.ras <- spatial_sync_raster(env.ras, template.ras, method="bilinear")
      env.ras <- mask(env.ras, template.ras)
      pos_df[,predictors$name[p]] <- extract(env.ras,pos_df[,c("long","lat")])
      cat("Predictor",predictors[p,"name"], "resampled and masked at resolution",resolution,"\n")

    } else {

      env.ras <- crop(env.ras, template.ras)
      stat <- predictors$resample[p]
      result <- special_resample(env.ras,template.ras, statistic=stat)
      result_df <- result$frame

      pos_df <- merge(pos_df, result$frame, all = F, all.x = T)
      colname <- paste(predictors$name[p],"_",stat,sep="")
      pos_df[,colname] <- pos_df$value
      colnums <- which(names(pos_df) %in% c("zone","cellID","value"))
      pos_df <- pos_df[, - colnums]

      cat("Predictor",predictors[p,"name"], "resampled and masked using function",stat,"\n")

    }
  }

  # filter to cells with all data
  pos_df <- na.omit(pos_df)

  # transform the response variable if needed (as set in the diversities data frame)
  if (diversities[j,"transform"] == "log") {
    pos_df$div <- log(pos_df$div)
  }

  # rescale the predictors - KEEP THE COLUMN NUMBERS CORRECT
  pos_df[,6:15] <- scale(pos_df[,6:15])

  # add a binary variable for membership of top class
  thresh <- quantile(pos_df$div,logistic_threshold)
  pos_df$div_binary <- rep(0,nrow(pos_df))
  pos_df$div_binary[pos_df$div >= thresh] <- 1

  # calculate the distance matrix for SARLM
  if (do_spatial_LM) {
    cat("\n\nPreparing neighbourhood weights for SARLM\n\n")
    coords<-as.matrix(cbind(pos_df$row, pos_df$col))
    cont.nb <- dnearneigh(coords,0,spatial_LM_radius,longlat=FALSE)
    weight_list <- nb2listw(cont.nb, glist=NULL, style="W", zero.policy=TRUE)
  }

  n <- nrow(pos_df)

  # automated predictor selection
  formula <- "div ~ bio1 + bio12 + bio17 + AET + Elev + Elev_sd + max_fpar_quantile + max_fpar_quantile:bio12"

  formula_LGM <- "div ~ bio1 + bio12 + bio17 + bio1DIFF + bio12PROP + bio17PROP + AET + Elev + Elev_sd + max_fpar_quantile + max_fpar_quantile:bio12"

  fullGLM      <- glm(data=pos_df, formula=formula)
  print(summary(fullGLM))
  r2 <- (fullGLM$null.deviance - fullGLM$deviance) / fullGLM$null.deviance

  fullGLM_LGM      <- glm(data=pos_df, formula=formula_LGM)
  print(summary(fullGLM_LGM))
  r2_LGM <- (fullGLM_LGM$null.deviance - fullGLM_LGM$deviance) / fullGLM_LGM$null.deviance

  cat("\nresolution",resolution,"\nr squared (present climate)",round(r2,4),"\nr squared (present  & LGM climate)", round(r2_LGM,4),"\n\n")

  chosen_GLM <- fullGLM_LGM ###############################################

  GLM_sum <- summary(chosen_GLM)
  #glm_one_res_results[,] <- glm_blank_results
  glm_one_res_results$resolution  <- resolution
  glm_one_res_results$formula     <- GLM_sum$bestmodel
  glm_one_res_results$r2          <- (chosen_GLM$null.deviance - chosen_GLM$deviance) / chosen_GLM$null.deviance
  glm_one_res_results$n           <- n

  # get model terms and co-efficents
  terms <- attributes(GLM_sum$terms)$term.labels
  coef  <- GLM_sum$coefficients
  for (pred in terms) {
    glm_one_res_results[1,pred] <- round(coef[pred,"Estimate"],4)
  }

  x <- glm_one_res_results

}

sfQuickStop()

print(scale_results)

file_path <- paste(results.dir,"GLM_by_scale_LGM.csv",sep="")
write.csv(scale_results, file_path, row.names=F)
rm(file_path)

# a plot
if (do_plots) {

  if (! exists("scale_results")) {
    scale_results       <- read.csv("C:/Users/u3579238/Work/AMT/Models/diversity/H_binoei/correlations/GLM_by_scale_LGM.csv")
    scale_results_NoLGM <- read.csv("C:/Users/u3579238/Work/AMT/Models/diversity/H_binoei/correlations/GLM_by_scale_no_LGM.csv")
    scale_results$r2_noLGM <- scale_results_NoLGM$r2
  }

  # plot model fit
  windows(10,12)
  plot(scale_results$resolution,scale_results$r2,xlab="resolution in degrees",cex.lab=1.2,ylab="GLM model r2", pch=20, ylim=c(0.1,1),log="x", col="red", type="p")
  points(scale_results$resolution,scale_results$r2_noLGM, pch=20, col="black")
  lines(lowess(scale_results$resolution, scale_results$r2, f=0.2),xlab="resolution in degrees",cex.lab=1.2, ylab="GLM model r2", ylim=c(0.1,1),log="x", col="red", lwd=1)
  lines(lowess(scale_results$resolution, scale_results$r2_noLGM, f=0.2), col="black", lwd=1)

  text(0.02,0.76,"black: model with present climate",col="black", pos=4)
  text(0.02,0.8,"red: model with present & LGM climate",col="red", pos=4)

  # now plot some parameters
  windows(10,12)
  plot(lowess(scale_results$resolution,scale_results$bio1,f=0.15),xlab="resolution in degrees",cex.lab=1.2,ylab="scaled parameter value", pch=20,log="x", col="red", type="l", ylim=c(-1.2,1.2), lwd=2)
  lines(lowess(scale_results$resolution,scale_results$bio1DIFF, f=0.15), col="red", lty=3, lwd=2)

  lines(lowess(scale_results$resolution,scale_results$bio12, f=0.15), col="blue", lwd=2)
  lines(lowess(scale_results$resolution,scale_results$bio12PROP, f=0.15), col="blue", lty=3, lwd=2)

  lines(lowess(scale_results$resolution,scale_results$bio17, f=0.15), col="green", lwd=2)
  lines(lowess(scale_results$resolution,scale_results$bio17PROP, f=0.15), col="green", lty=3, lwd=2)

  lines(lowess(scale_results$resolution,scale_results$AET, f=0.15), col="brown", lwd=2)

  abline(h = 0, col="black", lwd=2)

  text(0.015,1.1,"red: temperature",col="red", pos=4, cex=0.8)
  text(0.015,1.02,"red dotted: temperature (current - LGM)",col="red", pos=4, cex=0.8)
  text(0.015,0.94,"blue: annual precipitation",col="blue", pos=4, cex=0.8)
  text(0.015,0.86,"blue dotted: annual precipitation (current / LGM)",col="blue", pos=4, cex=0.8)
  text(0.015,0.78,"green: dry qtr precipitation",col="green", pos=4, cex=0.8)
  text(0.015,0.70,"green dotted: dry qtr precipitation (current / LGM)",col="green", pos=4, cex=0.8)
  text(0.015,0.62,"brown: AET",col="brown", pos=4, cex=0.8)

}