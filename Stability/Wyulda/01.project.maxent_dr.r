#drafted by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
#GNU General Public License .. feel free to use / distribute ... no warranties

library(raster)

################################################################################
#define directories
work.dir =        'C:/Users/u3579238/GISData/Helping/Sally/Wyulda_stability/maxent.output_2_5deg_topo_bathym/'; setwd(work.dir)
mxe.dir =    'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/mxe/'
maxent.jar = 'C:/Users/u3579238/Work/Refugia/Stability/maxent.jar'

# if clipping model outputs to a mask
mask.file =  'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/wyulda_buffer_2_5deg.asc'
clip_to_mask = TRUE

if (clip_to_mask) {
  mask.ras <- raster(mask.file)
}

species_name <- "Wyulda"

use_hpc = FALSE # If true, generic.skeletonrate hpc shell scripts to run models. If false, run maxent locally, one model at a time (DR 30 Aug 2012)

if (use_hpc) {
  tmp.sh.dir = 'C:/Users/Dan/Desktop/Dan_Craig_data/RF_Aug2012/tmp.sh.dir/'; dir.create(tmp.sh.dir); system(paste('rm -rf ',tmp.sh.dir,'/*',sep=''))
  setwd(tmp.sh.dir)
}

################################################################################
#list the projections, cycle thorugh them and project the models onto them
proj.list = list.files(mxe.dir) #list the projections

model_count = 0

#cycle through the projections
for (tproj in proj.list) {
  
  if (use_hpc) {
    
    ##create the sh file
    zz = file(paste(tproj,'.sh',sep=''),'w')
    cat('##################################\n',file=zz)
    cat('#!/bin/sh\n',file=zz)
    cat('cd $PBS_O_WORKDIR\n',file=zz)
    cat('java -cp ',maxent.jar,'density.Project ',work.dir,species_name,'.lambdas ',mxe.dir,tproj,' ',work.dir,tproj,'.asc fadebyclamping nowriteclampgrid\n',sep="",file=zz)
    cat('cd ',work.dir,'\n',sep='',file=zz)
    cat('gzip ',tproj,'.asc\n',sep='',file=zz)
    cat('##################################\n',file=zz)
    close(zz) 
    #submit the script
    system(paste('qsub -m n ',tproj,'.sh',sep=''))
    
  } else {
    model_count = model_count+1
    
    cat("\nAbout to project model for year", tproj,"\n")
    
    #Original model projection
    maxent_call = paste('java -mx1024m -cp ',maxent.jar,' density.Project ',work.dir,species_name,'.lambdas ',mxe.dir,tproj,' ',work.dir,tproj,'.asc fadebyclamping nowriteclampgrid',sep="")
    
    #Modified model projection to test a model fitted with a restricted set of predictors
    #maxent_call = paste('java -mx1024m -cp ',maxent.jar,' density.Project ',work.dir,species_name,'.lambdas ',mxe.dir,tproj,' ',work.dir,"/maxent.output1/",tproj,'.asc fadebyclamping nowriteclampgrid',sep="")
    
    cat(maxent_call,"\n")
    
    system(maxent_call, wait=TRUE) #run a model
    # NOTE 'wait=FALSE' in preceding line sends all jobs to the processors simulateously - uses full power of computer 
    # but could freeze other processes on fill memory (worked fine on my PC though)
  }
  
  if (clip_to_mask) {
    maxent_result.file = paste(work.dir,tproj,".asc",sep="")
    maxent_result.ras <- raster(maxent_result.file)
    maxent_result.ras <- mask(maxent_result.ras,mask.ras)
    writeRaster(maxent_result.ras,maxent_result.file,format="ascii",overwrite=TRUE)
    cat(maxent_result.file," clipped to mask\n")
  }
  
}

