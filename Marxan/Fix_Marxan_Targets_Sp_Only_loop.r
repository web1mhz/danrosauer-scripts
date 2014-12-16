rm(list=ls())

input_dir        <- "/home2/danr/marxan_mammals/input_100"
output_dir       <- "/home2/danr/marxan_mammals/input_100" # the output directory sits within this
dir_pattern      <- "run_"
input_filename   <- "spec_ph_25pc_cap.dat"
output_filename  <- 'spec_sp_25pc_cap.dat'

spf            <- 16.14367 # mean length of terminal branches across 100 trees

dirs = list.files(input_dir,pattern=dir_pattern, include.dirs=TRUE, recursive=FALSE)

# remove gz files from the vector of directories
with_gz <- grep("gz",dirs)

if (length(with_gz) > 0) {dirs <- dirs[- with_gz]}

for (dir in dirs) {

  filename_read <- paste(input_dir,dir,input_filename,sep="/")
  spec <- read.csv(filename_read)
  #names(spec) <- c("id","target","spf","name")
  terminals   <- which(substr(spec$name,1,4) != "node")

  #spec$target[-terminals] <- 0  # commented out to leave the targets unchanged, so performance against them can be easily assessed
  spec$spf[-terminals] <- 0

  # set the SPF for species
  spec$spf[terminals] <- spf
  
  #spec <- spec[,1:4] # omit the range column

  # write the version with new targets
  filename_new <- paste(output_dir,dir,output_filename,sep="/")
  write.csv(spec,filename_new,row.names=F)

  cat("\nDone for", dir, "\n")
}