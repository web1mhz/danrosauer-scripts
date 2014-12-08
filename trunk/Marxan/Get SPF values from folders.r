rm(list=ls())

base_dir       <- "/home2/danr/marxan_mammals/ph_25cap/" # the output directory sits within this
dir_pattern    <- "run_"
output_dir     <- base_dir
file_to_get   <- "spec_ph_25pc_cap.dat"

dirs = list.files(base_dir,pattern=dir_pattern, include.dirs=TRUE, recursive=FALSE)

result <- vector("numeric")
i <- 1

for (dir in dirs) {

  dir_num <- strsplit(dir,dir_pattern)[[1]][2]
  dir_num  <- as.numeric(dir_num)
  if (!is.na(dir_num)) {
    full_path = paste(base_dir,dir,"/",file_to_get,sep="")

    this.result <- read.csv(full_path)

    # drop rows for internal nodes
    this.result <- this.result[which(substr(this.result$name,1,4) != "node"),]
    result[i] <- mean(this.result$spf)

    cat("Mean spf for", dir, "is" ,result[i], "\n")

    i <- i+1
  }
}

cat("Overall mean spf from", i-1, "trees =", mean(result),"\n")