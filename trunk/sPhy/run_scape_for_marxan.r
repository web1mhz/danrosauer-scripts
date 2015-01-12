#May 2014. Matt Helmus, Caroline Tucker, Arne Mooers, Lanna Jin
##running "SCAPE" modules: analysis of phylogenetic metrics using simulated landscapes and trees.

##The simplified version omits the perturbations and missing information modules. Will add to later versions...

## Links to generate files for PD conservation scenarios in Marxan - Dan Rosauer

rm(list=ls())

##Load Libraries
require(plyr)
require(apTreeshape)
require(phytools)
require(ape)
require(geiger)
require(phangorn)
require(picante)
require(entropart)
require(phylobase)

#directory for source files
#setwd("~/Work/Software/danrosauer-scripts/sPhy")

##back end code for the scape script - generates landscapes using trees
source('~/Dropbox/sDiv S-Phy/groups/simulation/code/simulation/Working_Sims/gradient_cart_simulation_v3.1.R', chdir = TRUE)

##plotting and landscape perturbation functions
source('~/Dropbox/sDiv S-Phy/groups/simulation/code/simulation/Working_Sims/landscape_manipulation_plotting.R', chdir = TRUE)

##generate trees with a variation of shapes
source('~/Work/Software/danrosauer-scripts/sPhy/tree_generate_fct.R', chdir = TRUE)

##create Marxan input files for PD conservation scenario
source("~/Work/Software/danrosauer-scripts/sPhy/Prepare_Marxan_Inputs.r")


env_by_cell <- function(environmental_gradient,scape.size) {
  size <- scape.size + 1 # currently the actual number of rows / columns is scape.size + 1
  row_seq <- rep(environmentalgradient,times=size)
  col_seq <- rep(environmentalgradient,each=size)
  cell_env <- row_seq + col_seq
  return(cell_env)
}


# Number of sims to run
#  Each simulation has the same parameters, but unique:
#   * tree
#   * landscape
#   * marxan runs
number_of_sims <- 5

# tree parameters
tree_type = "rcoal"

# scape parameters
numsp            <-  64   #specify number of species, must be power of 2
scape.size       <-  39
g.center         <-  1
g.range          <-  1
g.repulse        <-  1
wd.all           <-  ((scape.size+1)^2)*0.02
signal.center    <-  TRUE
signal.range     <-  TRUE
same.range       <-  FALSE
repulse          <-  FALSE
center.scale     <-  1
range.scale      <-  1
repulse.scale    <-  1
site.stoch.scale <-  0
sd.center        <-  3
sd.range         <-  1
rho              <-  NULL
th               <-  5
maxabund         <-  10000
#subsample.percent <- 75
#detect.threshold  <- 5

# marxan parameters
base_dir          <- "C:/Users/u3579238/Work/Software/Marxan/scape_test/"
first_dir_num     <- 1
output_dir_prefix <- "run"
plot_each_run <- T    # produce a pdf with richness map and range-size histogram
target            <- 0.1  # proportion of the range of each branch to be protected
spf_multiplier    <- 10   # Marxan parameter for importance of meeting targets
                          # Marxan Species Protection Factor (SPF) = branch length x spf_multiplier

# set up a metadata frame
metadata <- data.frame(run_number = 0,
                       marxan_dir = "",
                       tree_type = "",
                       tree_colless = 0,
                       tree_gamma = 0,
                       numsp = 0,
                       scape.size = 0,
                       cell.count = 0,
                       g.center = 0,
                       g.range = 0,
                       g.repulse = 0,
                       wd.all = 0,
                       signal.center = FALSE,
                       signal.range =  FALSE,
                       same.range =  FALSE,
                       repulse =  FALSE,
                       center.scale = 0,
                       range.scale = 0,
                       repulse.scale = 0,
                       site.stoch.scale = 0,
                       sd.center = 0,
                       sd.range = 0,
                       th = 0,
                       species_with_range = 0,
                       occupied_cells = 0,
                       mean_range_size = 0,
                       mean_richness = 0,
                       quantile_10_species_ranges = 0,
                       date = "",
                       stringsAsFactors=FALSE)

# generate all trees and associated metadata
tree_gen<-tree_generate(numsp, number_of_sims, type="average") 	#generate trees
treeinfo<-tree_gen$treeinfo 	#save resulting tree info
treelist<-tree_gen$treelist		#save resulting trees
rm(tree_gen)

# loop to generate Marxan input files for each simulation
for (scenario_num in 1:number_of_sims) {

  input_tree <- treelist[[scenario_num]]

  scape_out<-scape(input_tree,
                   scape.size=scape.size,
                   g.center=g.center,
                   g.range=g.range,
                   g.repulse=g.repulse,
                   wd.all=wd.all,
                   signal.center=signal.center,
                   signal.range=signal.range,
                   same.range=same.range,
                   repulse=repulse,
                   center.scale = center.scale,
                   range.scale = range.scale,
                   repulse.scale = repulse.scale,
                   site.stoch.scale = site.stoch.scale,
                   sd.center=sd.center,
                   sd.range=sd.range,
                   rho=rho,
                   th=th)

  # note, rho currently has a value of NULL, and is not included in the metadata

  environmentalgradient<-seq(from=-1*(scape.size)/2,to=(scape.size)/2)
  cell_env <- env_by_cell(environmentalgradient,scape.size)

  ##use this output to create an abundance matrix and plots for the scape output.
  species<-colnames(scape_out$Y)
  abund_out<-landscape_abunds(scape_out,maxabund,species=species,figs=F)
  PA_mat<-abund_out$PA_mat

  #name the output folder and create if necessary
  if (number_of_sims > 1) {
    dir_number <- scenario_num + first_dir_num - 1  # allows for starting with a later number
    this_output_dir <- paste(base_dir,output_dir_prefix,"_",dir_number ,sep="")
  } else {
    this_output_dir <- paste(base_dir,output_dir_prefix,sep="")
  }
  if (! file.exists(this_output_dir)) {dir.create(this_output_dir)}

  input_tree  <- phylo4(input_tree)

  cat("\n\n**  Starting run", scenario_num, "**")
  cat("\nNumber of species:", numsp)
  cat("\nLandscape size:", scape.size + 1, "x", scape.size + 1)
  cat("\nCells:", (scape.size + 1)^2 )

  # remove zero range species
  species_ranges <- apply(PA_mat, MARGIN = 2, FUN = sum)
  zero_ranges <- which(species_ranges == 0)
  PA_mat <- PA_mat[,-zero_ranges]
  input_tree <- subset(input_tree, tips.exclude=zero_ranges)
  species <- species[-zero_ranges]
  species_ranges <- apply(PA_mat, MARGIN = 2, FUN = sum)

  # some summary numbers
  cell_count <- nrow(PA_mat)
  cells <- data.frame(cell = 1:cell_count,area = rep(1,cell_count))
  node_count  <- nrow(input_tree@edge)
  species_count <- length(species)

  richness       <- apply(PA_mat, 1, sum)
  occupied_cells <- which(richness > 0)
  occupied_cell_count <- length(occupied_cells)
  cat("\nOccupied_cells: ", occupied_cell_count, "\t\t", round(100 * occupied_cell_count / ((scape.size + 1)^2), 1), "%", sep="" )

  cat("\n\nSpecies with range:", species_count)

  # PDF of richness map and range size histogram for each simulation
  if (scenario_num == 1 | plot_each_run) {  # plot only first unless plot_each_run=TRUE

    pdf_filename <- paste(this_output_dir, "inputs.pdf", sep = "/")
    pdf(pdf_filename, 16, 8)

    par(mfrow=c(1,2))

    # species richness
    col_richnessmap <- heatcol <- colorRampPalette(c("white", "yellow", "red"))
    richness.matrix <- matrix(richness, nrow=scape.size+1, ncol=scape.size+1, byrow=TRUE)
    image(richness.matrix, col = terrain.colors(15), main="Species Richness")

    #species x area
    hist(species_ranges, 12, ylab="Number of species",xlab="Number of sites",main="Species Area Relationship",col="lightgrey")

    #tree
    input_tree_4d <- phylo4d(input_tree, tip.data=species_ranges)
    names(input_tree_4d@data) <- "Range size"
    treePlot(input_tree_4d)

    dev.off()
  }

  # expand the presence absence matrix to have a column for each internal node
  PA_mat_nodes <- matrix(nrow=cell_count,ncol=node_count)
  PA_mat_nodes[,1:species_count] <- PA_mat

  # give names to internal branches (added to species names)
  node_labels <- paste("node_", species_count+1:nNodes(input_tree), sep="")
  nodeLabels(input_tree) <- node_labels

  nodes <- getNode(input_tree,1:node_count)
  dimnames(PA_mat_nodes) <- list(NULL,names(nodes))
  BranchData  <- data.frame(BranchNames=names(nodes), BranchRange=as.numeric(rep(0,node_count)),stringsAsFactors=FALSE)
  BranchData$BranchRange[1:species_count] <- apply(PA_mat_nodes[,1:species_count],MARGIN=2,FUN=sum)

  # loop through the internal nodes to add an occurrence column for each to the presence/absence matric
  for (m in (species_count+1):node_count) {

    branch_name  <- names(descendants(input_tree, m, type="tips"))
    branch_occ    <- apply((PA_mat[,which(species %in% branch_name)]),MARGIN=1,FUN=max)
    PA_mat_nodes[,m] <- branch_occ
    BranchData$BranchRange[m] <- sum(branch_occ)

    node <- as.integer(nodes[m])

  }

  marxan_files_ok  <- MarxanInputs(
                      write.dir      = this_output_dir,
                      quad_list      = cells,
                      tree           = input_tree,
                      branch_data    = BranchData,
                      node_occ       = PA_mat_nodes,
                      target         = target,
                      spf_multiplier = spf_multiplier)

  metadata[scenario_num,] <-  c(run_number   =   scenario_num,
                                marxan_dir   =   this_output_dir,
                                tree_type    =   tree_type,
                                tree_colless = treeinfo[scenario_num,"Colless"],
                                tree_gamma   = treeinfo[scenario_num,"Gamma"],
                                numsp        =   numsp,
                                scape.size   =   scape.size,
                                cell_count   =   cell_count,
                                g.center,
                                g.range,
                                g.repulse,
                                wd.all,
                                signal.center,
                                signal.range,
                                same.range,
                                repulse,
                                center.scale,
                                range.scale,
                                repulse.scale,
                                site.stoch.scale,
                                sd.center,
                                sd.range,
                                th,
                                species_with_range = species_count,
                                occupied_cells = length(which(richness > 0)),
                                mean_range_size = mean(species_ranges),
                                mean_richness = mean(richness),
                                quantile_10_species_ranges = quantile(species_ranges, 0.1),
                                date  =  date()
                                )

  cat("\nScenario",scenario_num,"generated\n")


  gc() # garbage collect from memory

}

metadata_filename <- paste(base_dir,"metadata.csv",sep="")
write.csv(metadata, metadata_filename, row.names=F)

cat("\nMarxan files for ", number_of_sims, " simulations done!\n", date(),"\n", sep="")
