rm(list=ls())
library(gdata)

#source("/home2/danr/scripts/phylo_spatial/Prepare_Marxan_Inputs.r")

input_dir      <- "/home2/danr/marxan_mammals/ph_25cap/run_1"
base_dir       <- "/home2/danr/marxan_mammals/sp_25cap" # the output directory sits within this
file_to_change <- "spec_ph_25pc_cap.dat"
file_new_name  <- "spec_sp_25pc_cap.dat"
spf            <- 16.14367 # mean length of terminal branches across 100 trees

filename_read <- paste(input_dir,"/",file_to_change,sep="")
spec <- read.csv(filename_read)
names(spec) <- c("id","target","spf","name")

terminals <- which(substr(spec$name,1,4) != "node")

#spec$target[-terminals] <- 0  # commented out to leave the targets unchanged, so performance against them can be easily assessed
spec$spf[-terminals] <- 0

# set the SPF for species
spec$spf[terminals] <- spf
spec$name <- trim(spec$name)

# write the version with new targets
filename_new <- paste(base_dir,"/",file_new_name,sep="")
#spec <- spec[,1:4] # omit the range column
write.csv(spec,filename_new,row.names=F)

cat("\nDone for",base_dir,"\n")
