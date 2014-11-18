rm(list=ls())
library(gdata)

base_dir       <- "/home2/danr/marxan_mammals/ph_step_batch/" # the output directory sits within this
dir_pattern    <- "run__"
file_to_change   <- c("spec.dat")

dirs = list.files(base_dir,pattern=dir_pattern, include.dirs=TRUE, recursive=FALSE)

i <- 0
for (dir in dirs) {

  filename <- paste(base_dir,dir,"/",file_to_change,sep="")
    if (file.exists(filename)) {
      spec <- read.fwf(filename, widths=c(5,10,9,10))
      names(spec) <- c("id","target","spf","name")
      spec$target <- max(target * 5, 10)
      write.fwf(spec,filename)
    } else {
      cat("Did not find", filename,"\n")
    }
}
