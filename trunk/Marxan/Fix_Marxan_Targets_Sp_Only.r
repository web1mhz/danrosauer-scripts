rm(list=ls())
library(gdata)

source("/home2/danr/scripts/phylo_spatial/Prepare_Marxan_Inputs.r")

base_dir       <- "/home2/danr/marxan_mammals/sp_step/run_many_10pc" # the output directory sits within this
file_to_change <- "spec_10pc.dat"
file_new_name  <- "spec_10pc_sp_only.dat"

dir = base_dir
filename_read <- paste(dir,"/",file_to_change,sep="")
spec <- read.fwf(filename_read, widths=c(5,9,9,31),header=F,skip=1)
names(spec) <- c("id","target","spf","name")
#### here original range, than apply get_target function
spec$target[spec$id > 4777] <- 0
spec$spf[spec$id > 4777] <- 0
spec$spf[spec$id <= 4777] <- 1
spec$name <- trim(spec$name)

# write the version with new targets
filename_new <- paste(dir,"/",file_new_name,sep="")
spec <- spec[,1:4] # omit the range column 
write.csv(spec,filename_new,row.names=F)

cat("\nDone for",dir,"\n")
