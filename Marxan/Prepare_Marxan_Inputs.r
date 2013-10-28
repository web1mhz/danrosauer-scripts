# prepare Global Mammal Phylogenetic input files for Marxan
# Dan Rosauer - started September 2013

library(gdata)

MarxanInputs <- function(
                  write.dir,
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
  branch_targets  <- round(branch_ranges * target,4)  # for targets other than a flat %, set target via a function
  branch_info     <- as.data.frame(print(tree))
  branch_lengths  <- branch_info$edge.length
  branch_lengths[which(is.na(branch_lengths))] <- 0
  spec <- data.frame(cbind(branch_IDs,branch_targets,round(branch_lengths * spf_multiplier,4),branch_names))
  names(spec) <- c("species_id","target","spf","name")
  spec <- spec[-which(spec$spf==0),]

  pu_list  <- data.frame(sort(unique(node_occ[,2])))
  pu.count <- length(pu_list[,1])

  pu_list <- data.frame(cbind(rep(1:pu.count),pu_list))
  names(pu_list) <- c("pu_id","pu_name")

  pu <- data.frame(pu_list$pu_id)
  names(pu) <- "id"
  node_by_cell_SpID <- merge(spec, node_occ, by.x="species_id", by.y="BranchID")
  puvspr2 <- data.frame(node_by_cell_SpID$species_id,node_by_cell_SpID$QuadID,node_by_cell_SpID$Proportion)
  names(puvspr2) <- c("species","pu_name","amount")  # this data frame has the species number, but the PU name
  rm(node_by_cell_SpID)

  puvspr2 <- merge(pu_list, puvspr2, by="pu_name")
  puvspr2 <- puvspr2[,c(3,2,4)]
  names(puvspr2) <- c("species","pu","amount")
  
  # ensure correct column names for maxent format
  names(spec)[1] <- "id"
  names(pu) <- "id"

  try(dir.create(write.dir),silent = TRUE)

  setwd(write.dir)
  write.fwf(spec,"spec.dat")
  write.fwf(pu,"pu.dat")
  write.fwf(puvspr2,"puvspr2.dat")
  write.csv(pu_list,"pu_id_lookup.csv")
  
  return(1)
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
