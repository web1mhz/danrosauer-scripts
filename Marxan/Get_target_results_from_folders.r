rm(list=ls())

base_dir_ph       <- "/home2/danr/marxan_mammals/ph_25cap_10pc/" # the output directory sits within this
base_dir_sp       <- "/home2/danr/marxan_mammals/sp_25cap_10pc/" # the output directory sits within this
species_dir_ph    <- "/home2/danr/marxan_mammals/input_100/"
dir_pattern     <- "run_"
files_to_get    <- c("output_sum.txt","output_ssoln.txt")
file.pattern       <- 'output_mv0'
species_filename_ph   <- 'spec_ph_25pc_cap.dat'

target_set <- c("phylo","species")

for (TARGETS in target_set) {

  if (TARGETS == "species") {
    base_dir         <- base_dir_sp
  } else {
    base_dir         <- base_dir_ph
  }

  output_dir      <- paste(base_dir,"runs_combined/",sep="")

  species_filename <- species_filename_ph # use the spec file for the phylo version in all cases as it encodes the branch lengths to evaluate PD

  target_threshold   <- 0.95  # proportion of target considered 'meeting target' for threshold version

  dirs = list.files(base_dir,pattern=dir_pattern, include.dirs=TRUE, recursive=FALSE)

  i <- 0
  j <- 0

  for (dir in dirs) {

    dir_num <- strsplit(dir,dir_pattern)[[1]][2]
    dir_num  <- as.numeric(dir_num)
    if (!is.na(dir_num) & dir_num != 58) {  # temporary fix to crash in tree 58...
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
  	 solutions <- solutions[,c("Conservation.Feature","Target","MPM")]
          names(solutions)[3] <- paste("tree",dir_num,"_run",run,sep="")

        } else {
          # add a row or column for the next result
          solutions <- cbind(solutions,this.result[,"MPM"])
        }
        i <- i+1
        names(solutions)[i+2] <- paste("tree",dir_num,"_run",run,sep="")

        cat("Added tree",dir_num,"run",run,"\n")
      }

      # now load the species lookup file for that tree
      sp_file <- paste(species_dir_ph, dir, species_filename_ph, sep='/')
      sp_lookup <- read.csv(sp_file)

      if (j==0) {
        # create the initial data objects from the first run, as a weights (spf) df and a target df
        sp_weights <- sp_lookup
        sp_weights <- sp_weights[,c("id", "name", "spf")]
        names(sp_weights)[3] <- paste("tree",dir_num,sep="")

        sp_targets <- sp_lookup
        sp_targets <- sp_targets[,c("id", "name", "target")]
        names(sp_targets)[3] <- paste("tree",dir_num,sep="")

      } else {
        # add a row or column for the next result
        sp_weights <- cbind(sp_weights, sp_lookup$spf)
        names(sp_weights)[j+3] <- paste("tree",dir_num,sep="")

        sp_targets <- cbind(sp_targets, sp_lookup$target)
        names(sp_targets)[j+3] <- paste("tree",dir_num,sep="")
      }

      cat("Added target and spf for tree ",dir_num,"\n")
      j <- j+1
    }
  }

  # add the species / branch lookup
  solutions_sp <- merge(sp_lookup, solutions, by.x="id", by.y="Conservation.Feature")

  terminals   <- which(substr(sp_lookup$name,1,4) != "node")

  performance    <- solutions_sp[,-(1:5)]

  performance_PD     <- performance
  performance_PD_th  <- performance
  performance_sp     <- performance[terminals,]
  performance_sp_th  <- performance[terminals,]

  weight_sum_by_run  <- vector(mode="numeric", length=ncol(performance))

  for (j in 1:ncol(performance)) {

    treename     <- unlist(strsplit(names(performance_PD)[j], "_"))[1]
    weight_col   <- which(names(sp_weights) == treename)

    performance_PD[,j]    <- performance_PD[,j] * sp_weights[,weight_col]
    weight_sum_by_run[j]  <- sum(sp_weights[,weight_col])

    performance_PD_th[,j] <- (performance[,j] >= target_threshold) * sp_weights[,weight_col]

    performance_sp[,j]    <- performance_sp[,j]
    performance_sp_th[,j] <- (performance_sp_th[,j] >= target_threshold) * 1

  }

  # relative performance PD
  performance_PD_run         <- apply(performance_PD,2,sum)
  performance_PD_run_prop    <- performance_PD_run / weight_sum_by_run

  # threshold performance PD
  performance_PD_th_run      <- apply(performance_PD_th,2,sum)
  performance_PD_th_run_prop <- performance_PD_th_run / weight_sum_by_run

  # relative performance species
  performance_sp_run         <- apply(performance_sp,2,sum)
  performance_sp_run_prop    <- performance_sp_run / length(terminals)

  # threshold performance species
  performance_sp_th_run      <- apply(performance_sp_th,2,sum)
  performance_sp_th_run_prop <- performance_sp_th_run / length(terminals)

  if (TARGETS == "phylo") {
    results_ph <- data.frame(phylo_PD_prop   =  performance_PD_run_prop,
                          phylo_PD_thresh =  performance_PD_th_run_prop,
                          phylo_sp_prop   =  performance_sp_run_prop,
                          phylo_sp_thresh =  performance_sp_th_run_prop)
  } else {
    results_sp <- data.frame(species_PD_prop   =  performance_PD_run_prop,
                          species_PD_thresh =  performance_PD_th_run_prop,
                          species_sp_prop   =  performance_sp_run_prop,
                          species_sp_thresh =  performance_sp_th_run_prop)
  }

  setwd(output_dir)

  
  # draw histograms and scatterplots of performance separately for phylo and species optimizations
  
  pdf_name <- paste(TARGETS, "_performance_hist_",target_threshold,".pdf", sep="")
  pdf(pdf_name, 12, 12)

  par(mfrow=c(2,2), cex=1.2, cex.main=1)
  hist(performance_PD_run_prop, xlab="PD performance: branch length * proportion of target", main = "PD: Weighted by proportion of each target met")
  abline(v=mean(performance_PD_run_prop), col="blue", lwd=1.5)

  hist(performance_PD_th_run_prop, xlab="PD performance: branch length where target met", main=paste("PD: threshold ",target_threshold))
  abline(v=mean(performance_PD_th_run_prop), col="blue", lwd=1.5)

  hist(performance_sp_run_prop, xlab="Species performance: proportion of target", main = "Species: Weighted by proportion of each target met")
  abline(v=mean(performance_sp_run_prop), col="blue", lwd=1.5)

  hist(performance_sp_th_run_prop, xlab="Species performance: % where target met", main=paste("Species threshold: ",target_threshold))
  abline(v=mean(performance_sp_th_run_prop), col="blue", lwd=1.5)

  plot(performance_sp_run_prop,performance_PD_run_prop, xlab="Species performance", ylab="PD performance", pch=20, col="red", main="Species v PD captured (weighted)")
  plot(performance_sp_th_run_prop,performance_PD_th_run_prop, xlab="Species performance", ylab="PD performance", pch=20, col="red", main=paste("Species v PD captured (threshold ", target_threshold, ")", sep=""))

  dev.off()
  cat("\nNew plot saved for ", TARGETS, "\n")

}

# draw combined results
  
pdf_name <- paste("combined_performance_",target_threshold,".pdf", sep="")
pdf(pdf_name, 12, 12)
par(mfcol=c(2,2), cex=1.2, cex.main=1)

hist(results_ph$phylo_PD_prop, xlab="PD performance: branch length * proportion of target", main = "PD selected by PD: by proportion of each target") #, xlim=c(0.88,0.905))
abline(v=mean(results_ph$phylo_PD_prop), col="blue", lwd=2)

hist(results_sp$species_PD_prop, xlab="PD performance: branch length * proportion of target", main = "PD selected by species: by proportion of each target") #, xlim=c(0.88,0.905))
abline(v=mean(results_sp$species_PD_prop), col="blue", lwd=2)

hist(results_ph$phylo_PD_thresh, xlab="PD performance: branch length where target met", main=paste("PD selected by PD: threshold ", target_threshold)) #, xlim=c(0.76,0.815))
abline(v=mean(results_ph$phylo_PD_thresh), col="blue", lwd=2)

hist(results_sp$species_PD_thresh, xlab="Species performance: proportion of target", main = paste("PD selected by species: threshold ", target_threshold)) #, xlim=c(0.76,0.815))
abline(v=mean(results_sp$species_PD_thresh), col="blue", lwd=2)

dev.off()
cat("\nNew plot saved for species and phylo\n\n")

results <- cbind(results_ph, results_sp)


output_filename <- paste(output_dir,"marxan_sp_ph_performance.csv",sep="")
write.csv(results,output_filename,row.names=F)

