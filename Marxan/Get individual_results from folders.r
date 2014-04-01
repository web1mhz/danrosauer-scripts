rm(list=ls())

base_dir       <- "/home2/danr/marxan_mammals/ph_step_batch/" # the output directory sits within this
dir_pattern    <- "run__"
output_dir     <- paste(base_dir,"runs_combined/",sep="")
#files_to_get   <- c("pu.dat","input.dat")
files_to_get   <- c("output_sum.txt","output_ssoln.txt")
file.pattern       <- 'output_r'

runs_per_tree  <- 10

dirs = list.files(base_dir,pattern=dir_pattern, include.dirs=TRUE, recursive=FALSE)

i <- 0
for (dir in dirs) {

  dir_num <- strsplit(dir,dir_pattern)[[1]][2]
  dir_num  <- as.numeric(dir_num)
  if (!is.na(dir_num)) {
    full_path = paste(base_dir,dir,"/",sep="")
    output_files <- list.files(path = full_path, pattern = file.pattern, recursive = F,ignore.case = T, include.dirs = F)
    run <- 0
    
    for (file in output_files) {
      
      file <- paste(full_path,file,sep="")
      this.result <- read.csv(file)
      run <- run+1
  
      if (i==0) {
        # create the initial data objects from the first run, in two formats
        solutions <- this.result
        names(solutions)[2] <- paste("tree",dir_num,"_run",run,sep="")
        
        PU_names <- paste("PU",solutions$planning_unit,sep="")
        sol_alt_colnames <- c("tree","run",PU_names)
        sol_alt <- data.frame(t(this.result[,2]))
        sol_alt <- cbind(dir_num,run,sol_alt)
        
        names(sol_alt) <- sol_alt_colnames
        
      } else {
        # add a row or column for then next result
        solutions <- cbind(solutions,this.result[,2])
        sol_alt_new_row <- c(dir_num,run,this.result[,2])
        sol_alt   <- rbind(sol_alt,sol_alt_new_row)
      }
      i <- i+1
      names(solutions)[i+1] <- paste("tree",dir_num,"_run",run,sep="")
      
      cat("Added tree",dir_num,"run",run,"\n")
    }
  }
}

output_filename <- paste(output_dir,"all_solutions.csv",sep="")
write.csv(solutions,output_filename,row.names=F)

output_filename <- paste(output_dir,"all_solutions_alt_format.csv",sep="")
write.csv(sol_alt,output_filename,row.names=F)
