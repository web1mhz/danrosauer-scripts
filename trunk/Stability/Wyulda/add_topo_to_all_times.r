# reads in two specified .asc.gz layers and copies them as .asc.gz to the bioclim folder for each time period
# and donverts them to mxe in the mxe folder for each time period
# NOTE, this is set topo layers which do not change between time periods

rm(list=ls())

################################################################################
#define directories
source.dir = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000'
input.dir =  'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/mxe/000'
bioclim.dir = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/'
mxe.dir =    'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/mxe/'
maxent.jar = 'C:/Users/u3579238/Work/Refugia/Stability/maxent.jar'

################################################################################
#list the projections, cycle thorugh them 
proj.list = list.files(mxe.dir) #list the projections

setwd(source.dir)

#cycle through the projections
for (tproj in proj.list) {
  
  cat("\nCopying layers for year", tproj,"\n")
  folder <- paste(bioclim.dir,tproj,sep="")
  
  asc.file.1 <- paste(folder,"/slope_median_bathym.asc",sep="")
  asc.file.2 <- paste(folder,"/elev_range_bathym.asc",sep="")
  asc.file.old.1 <- paste(folder,"/slope_median.asc",sep="")
  asc.file.old.2 <- paste(folder,"/elev_range.asc",sep="")
  
  if (tproj != "000") {
    file.copy(from=paste(bioclim.dir,"000/slope_median_bathym.asc",sep=""), to=asc.file.1,overwrite=T)
    file.copy(from=paste(bioclim.dir,"000/slope_median_bathym.asc.gz",sep=""), to=paste(asc.file.1,".gz",sep=""),overwrite=T)
    file.copy(from=paste(bioclim.dir,"000/elev_range_bathym.asc",sep=""), to=asc.file.2,overwrite=T)
    file.copy(from=paste(bioclim.dir,"000/elev_range_bathym.asc.gz",sep=""), to=paste(asc.file.2,".gz",sep=""),overwrite=T)
  }
  
  # write to mxe
  maxent_call = paste('java -cp ',maxent.jar,' density.Convert ', bioclim.dir, tproj, '/ .asc ', mxe.dir, tproj, ' .mxe',sep="")
  cat(maxent_call,"\n")
 
  system(maxent_call, wait=TRUE, invisible=FALSE) #run a model
  
  if (tproj != "000") {
    file.remove(asc.file.old.1, asc.file.old.2)
    file.remove(paste(asc.file.1,".gz",sep=""),paste(asc.file.2,".gz",sep=""))
    mxe.old.1 <- paste(mxe.dir, tproj,'/elev_range.mxe',sep="")
    mxe.old.2 <- paste(mxe.dir, tproj,'/slope_median.mxe',sep="")    
    file.remove(mxe.old.1, mxe.old.2)
  }
}

