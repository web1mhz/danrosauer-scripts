rm(list=ls())
library(gdata)

base_dir       <- "/home2/danr/marxan_runs/mammals/batch_test/" # the output directory sits within this

start_dir <- getwd()

dirs = list.dirs(base_dir)

for (dir in dirs) {
  
  if ("sporder.dat" %in% list.files(dir)) { next }
  
  setwd(dir)
  puvspr2 <- read.fwf("puvspr2.dat")
  
  puvspr2 <- puvspr2[order(puvspr2$pu),]
  sporder <- puvspr2[order(puvspr2$species),]
  
  write.fwf(puvspr2,"puvspr2.dat")
  write.fwf(sporder,"sporder.dat")

  cat("\n",dir,date())
}
