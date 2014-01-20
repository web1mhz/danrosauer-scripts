rm(list=ls())

base_dir       <- "/home2/danr/marxan_mammals/batch_test/" # the output directory sits within this
dir_pattern    <- "marxan_batch_"
output_dir     <- paste(base_dir,"marxan_batch_combined/",sep="")
files_to_get <- c("pu.dat","input.dat")

dirs = list.files(base_dir,pattern=dir_pattern, include.dirs=TRUE, recursive=FALSE)

# combine results for output_ssoln.txt
fetch.filename="output_ssoln.txt"

i <- 0
for (dir in dirs) {

  dir_num <- strsplit(dir,"marxan_batch_")[[1]][2]
  dir_num  <- as.numeric(dir_num)
  if (!is.na(dir_num)) {
    filename <- paste(base_dir,dir,"/",fetch.filename,sep="")
      if (file.exists(filename)) {
        this.result <- read.csv(filename)
        if (i==0) {
          ssoln <- this.result
        } else {
          ssoln <- cbind(ssoln,this.result[,2])
        }
        i <- i+1
        names(ssoln)[i+1] <- paste("run",dir_num,sep="_")
      } else {
        cat("Did not find", filename,"\n")
      }
  } else {
    cat("Skipped",dir_num,"\n")
  }
}

# combine results for output_best.txt
fetch.filename="output_best.txt"

i <- 0
for (dir in dirs) {
  
  dir_num <- strsplit(dir,"marxan_batch_")[[1]][2]
  dir_num  <- as.numeric(dir_num)
  if (!is.na(dir_num)) {
    filename <- paste(base_dir,dir,"/",fetch.filename,sep="")
    if (file.exists(filename)) {
      this.result <- read.csv(filename)
      if (i==0) {
        best <- this.result
      } else {
        best <- cbind(best,this.result[,2])
      }
      i <- i+1
      names(best)[i+1] <- paste("run",dir_num,sep="_")
    } else {
      cat("Did not find", filename,"\n")
    }
  } else {
    cat("Skipped",dir_num,"\n")
  }
}


output_filename <- paste(output_dir,"ssoln.csv",sep="")
write.csv(ssoln,output_filename,row.names=F)

output_filename <- paste(output_dir,"best.csv",sep="")
write.csv(best,output_filename,row.names=F)