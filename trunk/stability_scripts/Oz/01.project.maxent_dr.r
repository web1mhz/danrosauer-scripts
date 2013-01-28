#drafted by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
#GNU General Public License .. feel free to use / distribute ... no warranties

################################################################################
#define directories
work.dir = 'C:/Users/u3579238/Work/Refugia/Stability/NE_NSW_RF/'; setwd(work.dir)
mxe.dir = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/mxe/'
maxent.jar = 'C:/Users/u3579238/Work/Refugia/Stability/maxent.jar'

use_hpc = FALSE # If true, generate hpc shell scripts to run models. If false, run maxent locally, one model at a time (DR 30 Aug 2012)

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
    cat('java -cp ',maxent.jar,'density.Project ',work.dir,'rf.lambdas ',mxe.dir,tproj,' ',work.dir,tproj,'.asc fadebyclamping nowriteclampgrid\n',sep="",file=zz)
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
    #maxent_call = paste('java -mx1024m -cp ',maxent.jar,' density.Project ',work.dir,'maxent.output/rf.lambdas ',mxe.dir,tproj,' ',work.dir,tproj,'.asc fadebyclamping nowriteclampgrid',sep="")
    
    #Modified model projection to test a model fitted with a restricted set of predictors
    maxent_call = paste('java -mx1024m -cp ',maxent.jar,' density.Project ',work.dir,'maxent.output1/rf.lambdas ',mxe.dir,tproj,' ',work.dir,"/maxent.output1/",tproj,'.asc fadebyclamping nowriteclampgrid',sep="")
    
    cat(maxent_call,"\n")
    
    system(maxent_call) #run a model
  }
}

