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
  #phylogeny_filename      <- "FritzTree.rs200k.100trees.tre"
  phylogeny_filename     <- "FritzTree.rs200k.100trees_1st_tree.tre"
  SpecMaster_filename     <- "Mammals 2010 Range Tree Lookup Nov2011.csv"
  base_dir                <- "/home2/danr/marxan_mammals/ph_25cap/" # the output directory sits within this
  output_dir              <- "runTest"
  output_type             <- "Marxan"
  land_cap                <- TRUE

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
  SpecMaster <- subset(SpecMaster, Family == "DASYURIDAE")

  #subset the species list to those in the spatial data, right up front
  SpecList   <- unique(Occ$SpecID)
  SpecMaster <- subset(SpecMaster, IUCNSpeciesID_Fritz %in% SpecList)

  #subset the species list terrestrial species
  SpecMaster <- subset(SpecMaster, Terrestrial == 1)

  #now subset the Occ data to the included species
  Occ <- Occ[Occ$SpecID %in% SpecMaster$IUCNSpeciesID_Fritz,]

  # combine duplicates due to multiple species joined to 1 tip on tree via SpecMaster
  Occ <- aggregate(Occ, by=list(Occ$QuadID, Occ$SpecID), FUN=sum)
  Occ <- Occ[,c("SpecID", "QuadID", "Proportion")]

  # ensure that the area occupied by a species in a cell is not greater than the cell area
  if (land_cap) {
    Occ_area <- merge(Occ,QuadInfo[,c("GRID360ID","LandProportion")], by.x="QuadID", by.y="GRID360ID", all=F)
    Occ_greater_land <- which(Occ_area$Proportion > Occ_area$LandProportion)
    Occ_area[Occ_greater_land,"Proportion"] <- Occ_area[Occ_greater_land,"LandProportion"]
    Occ <- Occ_area[,1:3]
    rm(Occ_area, Occ_greater_land)
  }


  #now subset the Occ data to the included species
  Occ <- Occ[Occ$SpecID %in% SpecMaster$IUCNSpeciesID_Fritz,]

  SpecMaster <- SpecMaster[,c(4,5)]
  names(SpecMaster) <- c("SpecID","taxon_name_tree")

  #number_of_cores <- 24
  #let doMPI and doMC work out the number of cores

  parameters1 <-  list(Trees_in = phy.orig,
                       QuadInfo_in = QuadInfo,
                      SpecMaster_in = SpecMaster,
                      Occ_in = Occ,
                      base_dir,
                      output_dir,
                      parallel_cores = 1,
                      multi_tree = TRUE,
                      iterations = 1,
                      interval = 1,
                      first_tree = 1,  # default is 1,
                      output_type,
			                first_dir_num = 1,
                      feedback = "Export node ranges")

  # run the PhyloSpatial function - observed
  results <- do.call('ExtractBranchRanges', parameters1)

  cat ("\n\n Started at:",StartTime,"\n")
  cat ("\n\n Finished at:",date(),"\n")
