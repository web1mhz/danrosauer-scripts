library(phylobase)

PhyloSpatialMulti = function(Trees_in,
                              QuadInfo_in, 
                              SpecMaster_in, 
                              Occ_in, 
                              base_dir,
                              output_dir,
                              parallel_cores = 1,
                              multi_tree = FALSE,
                              iterations = 1000,
                              interval = 1,
                              first_tree = 1,
                              output_Type = "generic",  #generic or Marxan
                              feedback = "")      #optional
{

  
   ##################
   #   PARAMETERS   #
   ##################
   
#  Tree_in is an ape phylo object
#  QuadInfo_in is a table which lists all quadrats (grid cells) by a column QuadID
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
  library(ape); library(sp);
    
  StartTime <- date()
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
  cat ("\n  Occurrence records loaded: ", length(Occ_in[,1]))
  Occ_in <- Occ_in[Occ_in$SpecID %in% SpecIDinTree,]
  cat ("\n  Occurrence records linked to tree were retained: ", length(Occ_in[,1]))  
    
  SpecList <- unique(Occ_in$SpecID)    # this way follows list and sequence as actually observed in GM360Occ
  TotalSpecCount <- length(SpecList); TotalSpecCount
  
  cat("\n\nStarting species calculations\n")
   ###############################
   # Species calculations
  
#    # range size
#   #QuadperSpec <- aggregate(Occ_in$SpecID, list(SpecID = Occ_in$SpecID), FUN = length)
#   QuadperSpec <- aggregate(Occ_in$SpecID, list(SpecID = Occ_in$SpecID), FUN = length)
#   QuadperSpec[,3] <- 1/QuadperSpec[2]
#   names(QuadperSpec) <- c("SpecID", "Range", "RangeWeight")
#   SpecAll_orig <- merge(SpecMaster, QuadperSpec, by.x="SpecID", by.y="SpecID", all.x=T)  
#   rm(QuadperSpec)
 
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
  block_size <- parallel_cores # this could be a number, or set to match the number of parallel cores
  
  loop_blocks <- seq(1,run_count,block_size)

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

      ###########################
      ### PhyloEndemism stuff ###
      ###########################
      phy4 <- phylo4(phy)
      CountNodes <- nrow(phy4@edge)
      
      # give names to internal branches (added to species names)
      node_labels <- paste("node_", 1:nNodes(phy4), sep="")
      nodeLabels(phy4) <- node_labels

      #################################################################
      cat("\nCalculating branch ranges for iteration",j,"\n")
  
      nodes <- getNode(phy4,1:CountNodes)
      #BranchIds   <- as.integer(nodes)
      BranchDone  <- rep(FALSE,CountNodes)
      BranchRange <- rep(0,CountNodes)
      BranchQuads <- rep(list(0),CountNodes)
      BranchData  <- data.frame(cbind(names(nodes), BranchDone, BranchRange, BranchQuads))
      phy4d       <- phylo4d(phy4,all.data=BranchData,merge.data=TRUE)
      tdata(phy4d)[,1] <- NULL
      rm(BranchDone,BranchRange,BranchQuads,BranchData)  # cleaning up
 
      for (m in 1:CountNodes) {
# a possible time saver:
  # check if BranchDone is true for all direct children.  
  # if so, just aggregate BranchQuads for the two, rather than for all descendents
                
        BranchLatin <- names(descendants(phy4d, m, type="tips")) 
        BranchSpecID <- SpecMaster$SpecID[SpecMaster$taxon_name_tree %in% BranchLatin]  #SpecIDs of species descendant from this node
        BranchOcc <- Occ[Occ$SpecID %in% BranchSpecID,]  # Subset Occ list to those species
        BranchQuad <- aggregate(BranchOcc[,2:3],by=list(BranchOcc$QuadID),FUN=sum)
        BranchQuad <- BranchQuad[,-2]
        names(BranchQuad) <- c("QuadID","Proportion")
        
        # ensure that where species occur in the same cell, the total proportion of area does not sum to > 1
        # including proportions of cells, rather than just whole cells, but without info on the location of species within each cell,
        # species ranges within a cell are treated as non-overlapping, until the total proportion of the cell occupied
        # exceeds 1.
        BranchQuad$Proportion[which(BranchQuad$Proportion>1)] <- 1
        # interestingly, without this line, I think we would be calculating AED

        tdata(phy4d)$BranchQuads[m] <- list(BranchQuad)
        tdata(phy4d)$BranchRange[m] <- sum(BranchQuad$Proportion)
        tdata(phy4d)$BranchDone[m] <- TRUE
        
        if (m %% 10 == 0) {
          cat(m,"of",CountNodes,"done\n")
        }
      }

#       cat("\n")  
#       BranchInfoAll <- cbind(BranchList, BranchLengths,  BranchRanges)
#       BranchOccAll <-  data.frame(cbind(BranchNumAll,BranchOccAll))
#       
#       names(BranchOccAll) <- c("BranchID","QuadID","Proportion")
#       BranchInfoAll <- data.frame(BranchInfoAll)
#       names(BranchInfoAll) <- c("BranchID","BranchName","BranchLength","BranchRS")
#       summary(BranchOccAll)
#       summary(BranchInfoAll)
      
      #QuadAll #item to store from each parallel loop
    }
  
    #name the output folder
    if (run_count > 1) {
      this_output_dir <- paste(base_dir,output_dir,"_",j,sep="")
    } else {
      this_output_dir <- paste(base_dir,output_dir,sep="")
    }

    if (output_Type == "Marxan") {
      result <- MarxanInputs(write.dir = this_output_dir,
                            tree_data = phy4d,
                            target =          0.1,
                            spf_multiplier =  3)
    }
# 
#     cat("\nAbout to write results for iterations", block_min,"to",block_max,"\n")
#     cat("\nIterations run: ", 1 + block_max - block_min," Available: ",length(unique(quad_scores[,1])),"\n")    
#     # Now write the results to file
#     if (k == 1 & new_file == TRUE) {
#       write.table(quad_scores, output_filename, col.names = TRUE, append = FALSE, row.names = FALSE, sep = ",", quote = FALSE,  eol = "\r\n", na = " ", dec = ".")
#     } else {
#       write.table(quad_scores, output_filename, col.names = FALSE, append = TRUE, row.names = FALSE, sep = ",", quote = FALSE,  eol = "\r\n", na = " ", dec = ".")
#     }  

  }

}

