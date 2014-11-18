rm(list=ls())
library(gdata)

base_dir       <- "/home2/danr/marxan_mammals/ph_step_batch/" # the output directory sits within this
dir_pattern    <- "run__"
file_to_change   <- c("sporder.dat")

dirs = list.files(base_dir,pattern=dir_pattern, include.dirs=TRUE, recursive=FALSE)

i <- 0
for (dir in dirs) {

  filename <- paste(base_dir,dir,"/",file_to_change,sep="")
    if (file.exists(filename)) {
      sporder <- read.fwf(filename, widths=c(5,6,7),skip=1)
      names(sporder) <- c("species","pu","amount")
      sporder <- sporder[order(sporder$species),]
      write.fwf(sporder,filename)
      cat(filename,"done\n")
    } else {
      cat("Did not find", filename,"\n")
    }
}
