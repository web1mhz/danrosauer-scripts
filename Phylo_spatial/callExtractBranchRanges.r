#! /usr/local/lib64/R/bin/Rscript --vanilla

  rm(list=ls())

  StartTime <- date()
  cat("\n\n************************************************\n")
  cat("Starting callExtractBranchRanges at: ", StartTime)
  cat  ("\n************************************************\n")
  
  source("/home2/danr/scripts/phylo_spatial/ExtractBranchRanges.r")
  source("/home2/danr/scripts/phylo_spatial/Prepare_Marxan_Inputs.r")
  setwd("/home2/danr/marxan_mammals/mammal_data")

  library(ape)
  
  cat("\n\nReading spatial files\n")

  # Files
  Occ_filename            <- "Occ360x114Mammals_2010_revised_Nov11.csv"
  #ClimLand_filename       <- "ClimLand_360x114_islands.csv"
  Global360x114_filename        <- "Global360x114.csv"
  phylogeny_filename      <- "FritzTree.rs200k.100trees.tre"
  #phylogeny_filename     <- "FritzTree.rs200k.100trees_1st_tree.tre"  
  SpecMaster_filename     <- "Mammals 2010 Range Tree Lookup Nov2011.csv"
  base_dir                <- "/home2/danr/marxan_mammals/ph_step_batch/" # the output directory sits within this
  output_dir              <- "run_"
  output_type             <- "Marxan"
  
  # Read and prepare occurrence file
  Occ <- read.csv(Occ_filename)
  Occ <- Occ[,c(5,2,10)]
  names(Occ) <- c("SpecID", "QuadID","Proportion")
  
  # Read alternative climland (Walter's revised) to get set of Quads to include  (with islands etc)
  QuadInfo  <- read.csv(Global360x114_filename)
  QuadInfo  <- subset(QuadInfo, select=c("GRID360ID", "X_COORD", "Y_COORD", "HBWID", "ROW", "COL", "AREA", "PERIMETER", "PROP0_0062", "LandProportion"))
  QuadInfo  <- subset (QuadInfo,  LandProportion >= 0.0002)
  Occ <- Occ[Occ$QuadID %in% QuadInfo[,1],]

  cat("\nReading phylogeny files\n")
  phy.orig <- read.nexus(file = phylogeny_filename)
  cat("\nTree loading finished\n")
  SpecMaster <- read.csv(SpecMaster_filename)
  cat("\n",length(SpecMaster[,1]),"\n")
 
  #TO SUBSET BY TAXA
  #SpecMaster <- subset(SpecMaster, Order == "DIPROTODONTIA")  

  #subset the species list to those in the spatial data, right up front
  SpecList   <- unique(Occ$SpecID)
  SpecMaster <- subset(SpecMaster, IUCNSpeciesID_Fritz %in% SpecList)

  #subset the species list terrestrial species
  SpecMaster <- subset(SpecMaster, Terrestrial == 1)
  
  #now subset the Occ data to the included species
  Occ <- Occ[Occ$SpecID %in% SpecMaster$IUCNSpeciesID_Fritz,]
  
  SpecMaster <- SpecMaster[,c(4,5)]
  names(SpecMaster) <- c("SpecID","taxon_name_tree")
  
  #phy <- phy.orig[[1]] # if several trees
  phy <- phy.orig
  
  #number_of_cores <- 24
  #let doMPI and doMC work out the number of cores
  
  parameters1 <-  list(Trees_in = phy,
                       QuadInfo_in = QuadInfo,
                      SpecMaster_in = SpecMaster, 
                      Occ_in = Occ, 
                      base_dir,
                      output_dir,
                      parallel_cores = 1,
                      multi_tree = TRUE,
                      iterations = 3,
                      interval = 1,
                      first_tree = 1,  # default is 1,
                      output_type,
			                first_dir_num = 1,
                      feedback = "Export node ranges")
#                     
#   if (TRUE) {
#     library(Rmpi)
#     if (mpi.comm.size(0) > 1) {
#       cat('must use mpirun -n 1\n')
#       mpi.quit()
#     }
#     
#     suppressMessages(library(doMPI))
#     cl <- startMPIcluster()
#     registerDoMPI(cl)
    
#     # run the PhyloSpatial function - spatial null model
#     if (exists("parameters2")) {
#       results <- do.call('PhyloSpatialMulti_sp_null_HPC', parameters2)  
#     }      
    # run the PhyloSpatial function - observed
    results <- do.call('ExtractBranchRanges', parameters1) 
    
    cat ("\n\n Started at:",StartTime,"\n")
#     cat ("\n\n Finished at:",date(),"\n")       
#     
#     closeCluster(cl)
# #     mpi.quit()
#   } else {
#     #suppressMessages(library(doMC))
#     
#     # for HPC clusters, let number of cores be set automatically
#     #registerDoMC(number_of_cores)
#     #registerDoMC(10)    
#     
#     # run the PhyloSpatial function
#     results <- do.call('PhyloSpatialMulti_sp_null_HPC', parameters2)  
#     results <- do.call('PhyloSpatialMulti', parameters1) 
#  }
