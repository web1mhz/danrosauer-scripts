rm(list=ls())

#base_dir       <- "/home2/danr/marxan_mammals/ph_25cap_10pc/" # the output directory sits within this
base_dir       <- "/home2/danr/marxan_mammals/sp_25cap_10pc/" # the output directory sits within this

dir_pattern    <- "run_"
output_dir     <- paste(base_dir,"runs_combined/",sep="")
#files_to_get   <- c("pu.dat","input.dat")
files_to_get   <- c("output_sum.txt","output_ssoln.txt")

runs_per_tree  <- 10

dirs = list.files(base_dir,pattern=dir_pattern, include.dirs=TRUE, recursive=FALSE)

# remove gz files from the vector of directories
with_gz <- grep("gz",dirs)
if (length(with_gz) > 0) {dirs <- dirs[- with_gz]}

# combine results for output_ssoln.txt
fetch.filename="output_ssoln.txt"

i <- 0
for (dir in dirs) {

  dir_num <- strsplit(dir,dir_pattern)[[1]][2]
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
    cat("Skipped",dir,"\n")
  }
}

# combine results for output_best.txt
fetch.filename="output_best.txt"

i <- 0
for (dir in dirs) {
  
  dir_num <- strsplit(dir,dir_pattern)[[1]][2]
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
    cat("Skipped",dir,"\n")
  }
}

output_filename <- paste(output_dir,"ssoln.csv",sep="")
write.csv(ssoln,output_filename,row.names=F)

output_filename <- paste(output_dir,"best.csv",sep="")
write.csv(best,output_filename,row.names=F)

# now combine to give a single value per cell
combined           <- data.frame(cbind(best[,1],apply(best[,-1],MARGIN=1,FUN=sum)),row.names=NULL)
names(combined)    <- c("planning_unit","best_sum")
combined$best_prop <- apply(best[,-1],MARGIN=1,FUN=mean)
combined$all_sum   <- apply(ssoln[,-1],MARGIN=1,FUN=sum)
combined$all_prop  <- apply(ssoln[,-1],MARGIN=1,FUN=mean)/runs_per_tree
combined$all_sd    <- apply(ssoln[,-1],MARGIN=1,FUN=sd)

output_filename <- paste(output_dir,"combined.csv",sep="")
write.csv(combined,output_filename,row.names=F)
