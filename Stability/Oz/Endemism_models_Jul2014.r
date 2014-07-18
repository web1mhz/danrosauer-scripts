
rm(list=ls())

library(SDMTools)
library(relaimpo)

# first a function - main script follows below
glm_process = function( result_frame,
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

base_path     <- 'C:/Users/u3579238/Work/Refugia/'
results.dir   <- paste(base_path,'Results/',sep='')

regions <- data.frame()
do_logistic <- TRUE
logistic_threshold <- 0.95  # for a logistic model to predict membership of top diversity

#whether to do spatial autocorrelation, and at what radius
do_spatial_LM <- FALSE
spatial_LM_radius <- 100

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

i <- 6
regions[i,"region"]        <- 'SEA'
regions[i,"veg_grid"]      <- paste(base_path,'Stability/NVIS/',regions[i,"region"],'_all_rainforest.asc',sep='')
regions[i,"model.dir"]     <- paste(base_path,'Stability/',regions[i,"region"],'_RF/',sep='')
regions[i,"stabil_static"] <- paste(regions[i,"model.dir"],'stability.topo.buf200km/static.sum.cost.asc',sep='')
regions[i,"stabil_10m"]    <- paste(regions[i,"model.dir"],'stability.topo.buf200km/shift.10.asc',sep='')
regions[i,"now_mod"]       <- paste(regions[i,"model.dir"],'maxent.output.topo.buf200km/000.asc',sep='')

diversities <- data.frame(taxon="lizardfrog",level="lineage", metric="endemism", grid="reptfrog_end_lin_25Sep_thresh_01.asc",stringsAsFactors = F)
j <- 2
diversities[j,1:4] <- c("lizardfrog","lineage","richness","reptfrog_rich_lin_25Sep_thresh_01.asc")
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
k      <- 0

# create a data frame to hold the model results
result_frame <- result_frame <- data.frame(region="",response_description="",predictor_description="",n=0,glm_r2=0,glm_aic=0,glm_delta_aic=0,sarlm_aic=0,sarlm_delta_aic=0,glm_logistic_r2=0,glm_logistic_aic=0,glm_logistic_delta_aic=0,sarlm_logistic_aic=0,sarlm_logistic_delta_aic=0,stringsAsFactors = F)

if (do_spatial_LM) {
  library(spdep)
}


#Loop through each region, and within it, each diversity metric
for (i in 1:nrow(regions)) {

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
    
    pos_rf = as.data.frame(which(is.finite(rf.asc),arr.ind=TRUE)) #get all points that have data
    pos_rf$rf = rf.asc[cbind(pos_rf$row,pos_rf$col)]            #append the vegetation data
    pos_rf$stabil_static = stabil_static.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$stabil_10m = stabil_10m.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$now_mod = now_mod.asc[cbind(pos_rf$row,pos_rf$col)] #append the stability data
    pos_rf$div = div_resample.asc[cbind(pos_rf$row,pos_rf$col)] #append the diversity data
    pos_rf = pos_rf[which(pos_rf$rf==1),]  # filter to rainforest areas
    pos_rf = pos_rf[which(pos_rf$div > 0),]  # filter to areas with a diversity score
    
    # add a binary variable for membership of top class
    thresh <- quantile(pos_rf$div,logistic_threshold)
    pos_rf$div_binary <- rep(0,nrow(pos_rf))
    pos_rf$div_binary[pos_rf$div >= thresh] <- 1

    # calculate the distance matrix for SARLM    
    if (do_spatial_LM) {
      cat("\n\nPreparing neighbourhood weights for SARLM\n\n")
      coords<-as.matrix(cbind(pos_rf$row, pos_rf$col))
      cont.nb <- dnearneigh(coords,0,spatial_LM_radius,longlat=FALSE)
      weight_list <- nb2listw(cont.nb, glist=NULL, style="W", zero.policy=TRUE)
    }
  
    n <- nrow(pos_rf)
    
    # fit models
    k <- k+1
    result_frame[k,"region"]  <- regions$region[i]
    result_frame$response_description[k] <- paste(diversities[j,1],diversities[j,2],diversities[j,3])
    result_frame$predictor_description[k] <- "Stability static"
    result_frame[k,"n"]       <- n
    result_frame <- glm_process(result_frame,
                                model_number=k,
                                predictor_text="stabil_static",
                                response_text = "div",
                                response_text_logistic="div_binary",
                                do_spatial_LM = do_spatial_LM,
                                weight_list = weight_list,
                                do_logistic=do_logistic,
                                the_data = pos_rf
                                )
    
    k <- k+1
    result_frame[k,"region"] <- regions$region[i]
    result_frame$response_description[k] <- paste(diversities[j,1],diversities[j,2],diversities[j,3])
    result_frame$predictor_description[k] <- "Stability 10m/yr"
    result_frame[k,"n"]       <- n  
    result_frame <- glm_process(result_frame,
                                model_number=k,
                                predictor_text="stabil_10m",
                                response_text = "div",
                                response_text_logistic="div_binary",
                                do_spatial_LM= do_spatial_LM,
                                do_logistic=do_logistic,
                                the_data = pos_rf
    )
    
    
    k <- k+1
    result_frame[k,"region"] <- regions$region[i]
    result_frame$response_description[k] <- paste(diversities[j,1],diversities[j,2],diversities[j,3])
    result_frame$predictor_description[k] <- "RF suitability now"
    result_frame[k,"n"]       <- n
    result_frame <- glm_process(result_frame,
                                model_number=k,
                                predictor_text="now_mod",
                                response_text = "div",
                                response_text_logistic="div_binary",
                                do_spatial_LM= do_spatial_LM,
                                do_logistic=do_logistic,
                                the_data = pos_rf
    )
  
  k <- k+1
  result_frame[k,"region"] <- regions$region[i]
  result_frame$response_description[k] <- paste(diversities[j,1],diversities[j,2],diversities[j,3])
  result_frame$predictor_description[k] <- "RF suitability now + Stability static"
  result_frame[k,"n"]       <- n
  result_frame <- glm_process(result_frame,
                              model_number=k,
                              predictor_text="now_mod + stabil_static",
                              response_text = "div",
                              response_text_logistic="div_binary",
                              do_spatial_LM= do_spatial_LM,
                              do_logistic=do_logistic,
                              the_data = pos_rf
  )
  
  k <- k+1
  result_frame[k,"region"] <- regions$region[i]
  result_frame$response_description[k] <- paste(diversities[j,1],diversities[j,2],diversities[j,3])
  result_frame$predictor_description[k] <- "RF suitability now + Stability 10m yr"
  result_frame[k,"n"]       <- n
  result_frame <- glm_process(result_frame,
                              model_number=k,
                              predictor_text="now_mod + stabil_10m",
                              response_text = "div",
                              response_text_logistic="div_binary",
                              do_spatial_LM= do_spatial_LM,
                              do_logistic=do_logistic,
                              the_data = pos_rf
  )
  #}
}

write.csv(output,paste(results.dir,"Stability_end_rich_cor_with_frogs_15Oct.csv",sep=""))

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

rm(rf.ras, rf.asc, div.ras, div.asc, div_resample.ras, stabil_static.asc, stabil_static.ras)

