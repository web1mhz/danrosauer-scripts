rm(list=ls())

base_dir       <- "/home2/danr/marxan_runs/mammals/batch_test/" # the output directory sits within this
source_dir    <- "/home2/danr/marxan_runs/mammals/batch_test/marxan_batch_1/"
files_to_copy <- c("pu.dat")

dirs = list.dirs(base_dir)

for (dir in dirs) {
  
  for (filename in files_to_copy) {
    file.copy(from=paste(source_dir,filename,sep=""),to=dir,overwrite=TRUE)
    cat(filename,"copied to",dir,"\n")
  }
}