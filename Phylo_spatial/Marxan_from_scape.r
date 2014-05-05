library(phylobase)

ExtractBranchRanges = function(Spmatrixlist,
                               treelist,
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
    
  StartTime <- date()
  cat("\n\n************************************************\n")
  cat("Starting PhyloSpatial at: ", StartTime)
  cat("\nUsing up to",parallel_cores,"cores for parallel stages")

  # get the trees from treelist
  for (tree in treelist) {
    tree <- as.phylo(tree[[1]])
  }


      if (output_Type == "Marxan") {
        result <- MarxanInputs(write.dir = this_output_dir,
                               quad_list = quad_list,
                               tree = phy4,
                               branch_data = BranchData,
                               node_occ = BranchOccAll,
                               target =          'func',
                               spf_multiplier =  3)
      }
      
      #QuadAll #item to store from each parallel loop

      try (rm(BranchData, BranchOccAll))
      gc() # garbage collect from memory
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