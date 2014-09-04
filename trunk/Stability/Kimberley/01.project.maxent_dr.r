#original version drafted by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
#GNU General Public License .. feel free to use / distribute ... no warranties

rm(list=ls())
library(raster)

################################################################################
#define directories
work.dir            <- 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/maxent.output_Carlia_johnstonei/'; setwd(work.dir)
mxe.dir             <- 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/mxe/'
maxent.jar          <- 'C:/Users/u3579238/Work/Refugia/Stability/maxent.jar'
mask_layer_name     <- 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/northern_mask.asc'  #only needed if

species_name        <- 'Carlia_johnstonei'

#clipping models to mask
clip_to_mask = TRUE

################################################################################
#list the projections, cycle thorugh them and project the models onto them
proj.list = list.files(mxe.dir) #list the projections

model_count = 0

if (clip_to_mask) {mask.ras <- raster(mask_layer_name)}

#cycle through the projections
for (tproj in proj.list) {

  model_count = model_count+1

  cat("\nAbout to project model for year", tproj,"\n")

  #Original model projection
  maxent_call = paste('java -mx1024m -cp ',maxent.jar,' density.Project ', work.dir, species_name, '.lambdas ',mxe.dir,tproj,' ',work.dir,tproj,'.asc fadebyclamping nowriteclampgrid',sep="")

  #Modified model projection to test a model fitted with a restricted set of predictors
  #maxent_call = paste('java -mx1024m -cp ',maxent.jar,' density.Project ',work.dir,'maxent.output1/C.johnstonei.lambdas ',mxe.dir,tproj,' ',work.dir,"/maxent.output1/",tproj,'.asc fadebyclamping nowriteclampgrid',sep="")

  cat(maxent_call,"\n")

  system(maxent_call) #run a model

  model.name <- paste(work.dir,tproj,".asc",sep="")
  mod.ras <- raster(model.name)

  if (clip_to_mask) {
    mod.ras <- crop(mod.ras,mask.ras)
    mod.ras <- mask(mod.ras,mask.ras)
    writeRaster(mod.ras,model.name,overwrite=TRUE)
  }

  plot(mod.ras,main=tproj)
}

