# prepare Global Mammal Phylogenetic input files for Marxan
# Dan Rosauer - started September 2013

library(gdata)

MarxanInputs <- function(
                  write.dir,
                  quad_list,
                  tree,
                  branch_data,
                  node_occ,
                  target,
                  spf_multiplier) {

  cat("\n\nGenerating marxan outputs for",write.dir,"\n\n")

  cons.feature.count <- nrow(tree@edge)
  branch_ranges   <- as.numeric(branch_data$BranchRange)
  branch_names    <- as.character(labels(tree))
  branch_IDs      <- as.numeric(getNode(tree))

  if (target > 1) {
    branch_targets <- round(get_target(branch_ranges,target),4)
  } else {
    branch_targets  <- round(branch_ranges * target,4)
  }

  branch_info     <- as.data.frame(print(tree))
  branch_lengths  <- branch_info$edge.length
  branch_lengths[which(is.na(branch_lengths))] <- 0
  spec <- data.frame(cbind(branch_IDs,branch_targets,round(branch_lengths * spf_multiplier,4),branch_names), stringsAsFactors=F)
  names(spec) <- c("species_id","target","spf","name")
  spec <- spec[-which(spec$spf==0),]

  # get the mean spf to define the spf for the species version
  species_count <- length(tipLabels(tree))
  mean_spf      <- mean(as.double(spec$spf[1:species_count]))

  pu_from_occ <- data.frame(unique(node_occ[,2]))
  pu_from_occ <- data.frame(QuadID = sort(pu_from_occ$QuadID))
  #names(pu_from_occ)  <- "QuadID"
  names(quad_list)[1] <- "QuadID"

  pu_list      <- merge(pu_from_occ,quad_list,by="QuadID", all=F)
  pu.count     <- nrow(pu_list)

  pu_list <- data.frame(cbind(rep(1:pu.count),pu_list))
  names(pu_list) <- c("pu_id","pu_name","LandProportion")

  pu <- data.frame(pu_list[,c(1,3)])
  names(pu) <- c("id","cost")
  node_by_cell_SpID <- merge(spec, node_occ, by.x="species_id", by.y="BranchID")

  puvspr2 <- node_by_cell_SpID[,c("species_id","QuadID","Proportion")]
  names(puvspr2) <- c("species","pu_name","amount")  # this data frame has the species number, but the PU name
  rm(node_by_cell_SpID)

  puvspr2 <- merge(pu_list, puvspr2, by="pu_name")
  puvspr2 <-  puvspr2[,c("species","pu_id","amount")]

  names(puvspr2) <- c("species","pu","amount")
  puvspr2$species <- as.numeric(puvspr2$species)

  # now order files to use Marxan quick preparation method
  puvspr2 <- puvspr2[order(puvspr2$pu),]
  sporder <- puvspr2[order(puvspr2$species),]

  # ensure correct column names for maxent format
  names(spec)[1] <- "id"

  try(dir.create(write.dir),silent = TRUE)

  cat("\nAbout to write files to ",write.dir,"\n")

  setwd(write.dir)
  write.csv(spec,"spec.dat", row.names=F)
  write.csv(pu,"pu.dat", row.names=F)
  write.csv(puvspr2,"puvspr2.dat", row.names=F)
  write.csv(sporder,"sporder.dat", row.names=F)
  write.csv(pu_list,"pu_id_lookup.csv", row.names=F)

  return(mean_tip_spf=mean_spf)
}

############################################

fixnames <- function (text_in) {
  #   if (grep("___",text_in)==1) {
  #     text_in <- paste("b",text_in,sep="")
  #   }
  text_out <- gsub('___',"",text_in)
  text_out <- gsub(" ","_",text_out)
  text_out <- paste("x_",text_out,sep="")
  return(text_out)
}

get_target <- function(range, option=1) {

# returns a target area, according to a function (this could be made more general with more parameters)
  # option 1: target = 10% of range
  # option 2: target = 50% if < 20 cells, 10 cells if > 20 cells
  # option 3: target = 25% if less than 100 cells, then a flat 25.  Minimum target = 1
  # option 3: target = 100% if < 10 cells, 10 cells if < 100, 10% if < 2000, then flat at 200.

  if (!(option %in% c(1,2,3,4))) {
    option <- 1
    cat("\nDefaulting to target option 1\n")
  }

  n <- length(range)
  target <- vector(mode="numeric",length=n)
  if (option ==1) {
    target <- range * 0.1
  } else if (option ==2) {
    for (i in 1:n) {
      target[i] <- min(range[i] * 0.5,10)
    }
  } else if (option == 3) {
    for (i in 1:n) {
      target[i] <- min(range[i] * 0.25,25)
      target[i] <- max(target[i], 1)
    }
    target[range < 1] <- range[range < 1]
  } else if (option ==4) {
    for (i in 1:n) {
      if (range[i] >= 100) {
        target[i] <- min(range[i],10)
      } else if (range[i] > 100) {
        target[i] <- min((range[i] / 10),200)
      }
    }
  }
  return(target)
}
