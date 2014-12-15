library(phylobase)
library(dplyr)

ExtractBranchRanges = function(Trees_in,
                              QuadInfo_in,
                              SpecMaster_in,
                              Occ_in,
                              base_dir,
                              output_dir,
                              parallel_cores = 1,
                              multi_tree = FALSE,
                              iterations = 100,
                              interval = 1,
                              first_tree = 1,
                              output_Type = "generic",  #generic or Marxan
  	  	                      first_dir_num = 1,
                              feedback = "")      #optional
{


   ##################
   #   PARAMETERS   #
   ##################

#  Tree_in is an ape phylo object
#  SpecMaster_in lists all species in the analysis with columns as follows:
#    taxon_name_tree gives taxon names which match the tree tip labels
#    SpecID gives a number (or name) to match values in Occ_in
#  Occ_in lists species occurence in each quadrat as follows:
#    SpecID matches a species in SpecMaster_in
#    QuadID matches a quadrat in QuadInfo_in
#  parallel_cores is an optional parameter to specify how many cores to use for
#     parallel processing (default = 4)
#  iterations specifies the number of times to repeat the randomization
#  first_tree is an optional parameter to all for restarting an analysis mid way
#     through the set of randomizations

  #library(doMC)
  library(ape);

  StartTime    <- date()
  StartSysTime <- Sys.time()

  cat("\n\n************************************************\n")
  cat("Starting PhyloSpatial at: ", StartTime)
  cat("\nUsing up to",parallel_cores,"cores for parallel stages")

  # add an optional feedback line to track progress
  if (feedback != "") {
    cat("\n",feedback, sep = "")
  }
  cat  ("\n************************************************\n")

  if (parallel_cores == 0) {parallel_cores <- 1}

  #registerDoMC(parallel_cores)

  ##-- input data --
  # check if Trees_in is a single tree or a set of trees
  if (class(Trees_in) == "multiPhylo") {
    phy <- Trees_in[[first_tree]] #use the first tree for name matching etc
  } else {
    phy <- Trees_in
  }
  SpecMaster <- SpecMaster_in

  ##-------- Tree calculations ------------------------------------------
  cat("\nDoing tree calculations\n")

  cat("\nStarting counts\n")
  # Counts
  cat("\n  Count of species or tips in the tree: ",length(phy$tip.label))
  cat("\n  Count of species in the spatial data: ",length(SpecMaster$taxon_name_tree)) # count of species in the data
  cat("\n  Count of species in spatial data that are in tree: ",length(intersect(phy$tip.label, as.character(SpecMaster$taxon_name_tree)))) # count of species in data that are in tree
  SpecIDinTree <- SpecMaster$SpecID [is.element(as.character(SpecMaster$taxon_name_tree), phy$tip.label)] # SpecIDs of species in data that are in tree
  cat("\n  Count of species in spatial data that are not in tree:",length(setdiff(as.character(SpecMaster$taxon_name_tree), phy$tip.label))) # count of species in data that are not in tree
  unmatched_tips <- setdiff(phy$tip.label, as.character(SpecMaster$taxon_name_tree))  # vector of species in tree that are not in data
  cat("\n  Count of species in tree that are not in spatial data",length(unmatched_tips),"\n")

  # trim the first tree to species that are in the spatial data, if required
  branch_count_orig <- length(phy$edge[,1])
  if (length(unmatched_tips) >= 1) {
    phy <- drop.tip(phy,unmatched_tips, trim.internal=TRUE)
    cat("\n  Tips and internal branches not represented in the spatial data have been removed.\n")
    cat("  Branches original: ", branch_count_orig,"\n")
    cat("  Branches remaining: ", length(phy$edge[,1]),"\n")
  } else {
    cat("\n  All", branch_count_orig ,"tree branches successfully matched to spatial data.\n")
  }

  # subset Occ to species that are in the tree
  cat ("\n  Occurrence records loaded: ", length(Occ_in[,1]),"\n")
  Occ_in <- Occ_in[Occ_in$SpecID %in% SpecIDinTree,]
  cat ("\n  Occurrence records linked to tree were retained: ", length(Occ_in[,1]),"\n")

  SpecList <- unique(Occ_in$SpecID)    # this way follows list and sequence as actually observed in GM360Occ
  TotalSpecCount <- length(SpecList); TotalSpecCount

  # tree preparation outside of the loop
  TipsCount <- length(phy$tip.label)
  reorder(phy, order = "cladewise")

  if (multi_tree) {
    trees_to_use <- seq(first_tree,first_tree + ((iterations-1)*interval),interval)
  } else {
    trees_to_use <- rep(first_tree, iterations)
  }

  # i is the tree number, j is the iteration number
  run_count <- length(trees_to_use)

  # block size defines how many analyses to do in the parallel loop before writing the results to disk
  # larger should be more efficient, but smaller is lees likely to blocw the memory
  block_size  <- parallel_cores # this could be a number, or set to match the number of parallel cores

  loop_blocks <- seq(1,run_count,block_size)

  spf.df      <- vector("double", run_count) # a vector to record the mean species SPF from each run

  #START LOOP THROUGH ALL TREES HERE, IN BLOCKS
  for (k in seq(along.with=loop_blocks)) {

    block_min <- loop_blocks[k]
    if (loop_blocks[k] == max(loop_blocks)) {
      block_max <- run_count
    } else {
      block_max <- loop_blocks[k+1]-1
    }

    #START PARALLEL LOOP THROUGH A BLOCK OF TREES HERE
    #quad_scores <- foreach (j = block_min:block_max, .combine = "rbind", .packages = "ape") %dopar%{
    for (j in block_min:block_max) {

      i <- trees_to_use[j]

      cat("\n********************************************************")
      cat("\n  Starting iteration", j, "using tree ", i)
      cat(" at: ",date())
      # an optional feedback line to track progress
      if (feedback != "") {
        cat("\n ",feedback,"\n")
      }
      cat("********************************************************\n")

      #use the same version of Occ and SpecAll for each iteration
      Occ <- Occ_in
      #SpecAll <- SpecAll_orig

      #trim unused tips and branches to species that are in the spatial data, if required
      if (i > first_tree) {
        phy <- Trees_in[[i]]
        branch_count_orig <- length(phy$edge[,1])
        if (length(unmatched_tips) >= 1) {
          phy <- drop.tip(phy,unmatched_tips, trim.internal=TRUE)
          cat("  Branches original: ", branch_count_orig,"\n")
          cat("  Branches after trim to match spatial data: ", length(phy$edge[,1]),"\n")
        } else {
          cat("\n  All", branch_count_orig ,"tree branches successfully matched to spatial data.\n")
        }
      }

      phy4 <- phylo4(phy)
      CountNodes <- nrow(phy4@edge)

      # give names to internal branches (added to species names)
      node_labels <- paste("node_", 1:nNodes(phy4), sep="")
      nodeLabels(phy4) <- node_labels

      #################################################################
      cat("\nCalculating branch ranges for iteration",j,"\n")
      nodes <- getNode(phy4,1:CountNodes)
      BranchDone  <- rep(FALSE,CountNodes)
      BranchRange <- as.numeric(rep(0,CountNodes))
      BranchData  <- data.frame(cbind(names(nodes), BranchDone, BranchRange),stringsAsFactors=FALSE)
      names(BranchData)[1] <- "BranchNames"
      rm(BranchDone,BranchRange)  # cleaning up

      tree_table <- as.data.frame(print(phy4))
      is.tip <- tree_table$node.type == "tip"

      for (m in 1:CountNodes) {

        if (is.tip[m]) {
          BranchLatin  <- tree_table$label[m]
          BranchSpecID <- SpecMaster$SpecID[SpecMaster$taxon_name_tree == BranchLatin]  #SpecIDs of species descendant from this node
          BranchOcc    <- Occ[Occ$SpecID %in% BranchSpecID,]  # Subset Occ list to those species
          BranchQuad   <- BranchOcc[,-2]

        } else {

          ######### TESTING
          children <-  tree_table$node[tree_table$ancestor==m]
          children_done <- BranchData$BranchDone[children]
	   if (all(children_done)) {
            childrenOcc  <- BranchOccAll[BranchOccAll$BranchID %in% children,]

            childrenOcc  <- summarise(group_by(childrenOcc[,2:3], QuadID), Proportion=sum(Proportion))

            BranchQuad   <- childrenOcc
            cat(".")  # remove once its working
	   } else {

          ######### TESTING

            BranchLatin  <- names(descendants(phy4, m, type="tips"))
            BranchSpecID <- SpecMaster$SpecID[SpecMaster$taxon_name_tree %in% BranchLatin]  #SpecIDs of species descendant from this node
            BranchOcc    <- Occ[Occ$SpecID %in% BranchSpecID,]  # Subset Occ list to those species
            #BranchQuad  <- aggregate(BranchOcc[,c(1,3)],by=list(BranchOcc$QuadID),FUN=sum)
            BranchQuad   <- summarise(group_by(BranchOcc[,c(1,3)], QuadID), Proportion=sum(Proportion))
            #BranchQuad   <- BranchQuad[,-2]
          }

          # ensure that where species occur in the same cell, the total proportion of area does not sum to > 1
          # including proportions of cells, rather than just whole cells, but without info on the location of species within each cell,
          # species ranges within a cell are treated as non-overlapping, until the total proportion of the cell occupied
          # exceeds 1.
          BranchQuad$Proportion[which(BranchQuad$Proportion>1)] <- 1

        }

        node <- as.integer(nodes[m])
        BranchQuad     <- cbind(rep(node,nrow(BranchQuad)),BranchQuad)
        names(BranchQuad) <- c("BranchID","QuadID","Proportion")

        # interestingly, without this line, I think we would be calculating AED

        if (m==1) {
          BranchOccAll <- BranchQuad
        } else {
          BranchOccAll <- rbind_list(BranchOccAll, BranchQuad)
          #BranchOccAll <- rbind(BranchOccAll, BranchQuad)
        }

        BranchData$BranchRange[m] <- sum(BranchQuad$Proportion)
        BranchData$BranchDone[m] <- TRUE

        # progress output - remove once working
        if (m %% 100 == 0) {
          cat("\n",m,"of",CountNodes,"done for tree\t",i)
	   minutes_elapsed <- as.double(round(difftime(Sys.time(), StartSysTime, units = "mins"), 2))
          cat("\n",tree_table$label[m],"\t",nrow(BranchQuad),"\t",nrow(BranchOccAll),"\t", ( j + first_dir_num - 1), "\t", minutes_elapsed, "minutes\n")
        }
      }

      #name the output folder
      if (run_count > 1) {
        this_output_dir <- paste(base_dir, output_dir, "_", i ,sep="")  # use the tree number as the directory number
      } else {
        this_output_dir <- paste(base_dir, output_dir, sep="")
      }

      quad_list <- QuadInfo_in[,c("GRID360ID", "LandProportion")]

      gc()

##############################################################
#### TEMPORARY FEEDBACK ####
cat("\nRows in BranchOccAll:",nrow(BranchOccAll),"\n")
cat("\nUnique rows in BranchOccAll:",nrow(unique(BranchOccAll)),"\n")
cat("\n",names(BranchOccAll)[1:2],"\n")
cat("\nUnique BranchOccAllin Occ cols 1-2:",nrow(unique(BranchOccAll[,1:2])),"\n")
#### TEMPORARY FEEDBACK ####

      if (output_Type == "Marxan") {
        result <- MarxanInputs(write.dir = this_output_dir,
                               quad_list = quad_list,
                               tree = phy4,
                               branch_data = BranchData,
                               node_occ = BranchOccAll,
                               target =          3, # 25% if < 100 cells, then a flat 25. Min target = 1
                               spf_multiplier =  3)
      }

      if (length(spf.df) > 1) {
        spf.df[j] <- result
      } else {
        spf.df    <- result
      }

      #QuadAll #item to store from each parallel loop

      try (rm(BranchData, BranchOccAll))
      gc() # garbage collect from memory
    }

  }

  # store the mean SPF values, to create a species Marxan run with SPF matching the mean value for all species
  write.csv(spf.df,"mean_species_spf.csv", row.names=F)

}