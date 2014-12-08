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
  library(plyr)
  
  cat("\n\nReading spatial files\n")

  # Files
  Occ_filename            <- "Occ360x114Mammals_2010_revised_Nov11.csv"
  Global360x114_filename        <- "Global360x114.csv"
  phylogeny_filename      <- "FritzTree.rs200k.100trees.tre"
 #phylogeny_filename      <- "FritzTree.rs200k.100trees_1st_tree.tre"
  SpecMaster_filename     <- "Mammals 2010 Range Tree Lookup Nov2011.csv"
  base_dir                <- "/home2/danr/marxan_mammals/ph_25cap/" # the output directory sits within this
  output_dir              <- "runTest"
  output_type             <- "Marxan"
  land_cap                <- TRUE
  first_tree		     <- 2
  number_of_trees         <- 1

  # Read and prepare occurrence file
  Occ <- read.csv(Occ_filename)
  #Occ <- Occ[,c(5,2,10)]
  Occ <- Occ[,c(8,2,10)]  # Using FritzSpeciesID as the SpecID, rather than IUCN_Taxon_ID_Fritz
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
  cat("\nRecords in species list:",length(SpecMaster[,1]),"\n")

  #TO SUBSET BY TAXA
  #SpecMaster <- subset(SpecMaster, Order == "DIPROTODONTIA")

  #subset the species list to those in the spatial data, right up front
  SpecList   <- unique(Occ$SpecID)
  SpecMaster <- subset(SpecMaster, Fritz_species_ID %in% SpecList)

  #subset the species list terrestrial species
  SpecMaster <- subset(SpecMaster, Terrestrial == 1)

  #now subset the Occ data to the included species
  #Occ <- Occ[Occ$SpecID %in% SpecMaster$IUCNSpeciesID_Fritz,]
  Occ <- Occ[Occ$SpecID %in% SpecMaster$Fritz_species_ID ,] # Using FritzSpeciesID as the SpecID, rather than IUCN_Taxon_ID_Fritz

##################################################################
#  # replace Occ$SpecID with Fritz_species_ID to ensure that anything linked to the same tip has the same specID
#  SpecIDTranslate <- SpecMaster[,c("IUCNSpeciesID_Fritz","Fritz_species_ID")]
#  Occ_merge <- merge(Occ, SpecIDTranslate, by.x="SpecID", by.y="IUCNSpeciesID_Fritz")
#  Occ_merge$SpecID <- Occ_merge$Fritz_species_ID
#  Occ <- Occ_merge[,1:3]
#  rm(Occ_merge)
##################################################################

  cat("\nStarting ddply for", nrow(Occ), "rows in Occ at ", date(),"\n")
  Occ <- ddply(Occ, .(QuadID, SpecID), summarize, Prop=sum(Proportion))
  names(Occ)[3] <- "Proportion"
  cat("\nFinished ddply with", nrow(Occ), "rows in Occ at ", date(),"\n")

  # ensure that the area occupied by a species in a cell is not greater than the cell area
  if (land_cap) {
    Occ_area <- merge(Occ,QuadInfo[,c("GRID360ID","LandProportion")], by.x="QuadID", by.y="GRID360ID", all=F)
    Occ_greater_land <- which(Occ_area$Proportion > Occ_area$LandProportion)
    Occ_area[Occ_greater_land,"Proportion"] <- Occ_area[Occ_greater_land,"LandProportion"]
    Occ <- Occ_area[,1:3]
    rm(Occ_area, Occ_greater_land)
  }

  SpecMaster <- unique(SpecMaster[,c(6,5)])
  names(SpecMaster) <- c("SpecID","taxon_name_tree")

  number_of_cores <- 1
  #let doMPI and doMC work out the number of cores

  parameters1 <-  list(Trees_in = phy.orig,
                       QuadInfo_in = QuadInfo,
                      SpecMaster_in = SpecMaster,
                      Occ_in = Occ,
                      base_dir,
                      output_dir,
                      parallel_cores = 1,
                      multi_tree = TRUE,
                      iterations = number_of_trees,
                      interval = 1,
                      first_tree = first_tree,  # default is 1,
                      output_type,
			 first_dir_num = 1,
                      feedback = "Export node ranges")

  # run the PhyloSpatial function - observed
  results <- do.call('ExtractBranchRanges', parameters1)

  cat ("\n\n Started at:",StartTime,"\n")
  cat ("\n\n Finished at:",date(),"\n")
