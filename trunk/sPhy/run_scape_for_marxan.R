#May 2014. Matt Helmus, Caroline Tucker, Arne Mooers, Lanna Jin
##running "SCAPE" modules: analysis of phylogenetic metrics using simulated landscapes and trees.

##The simplified version omits the perturbations and missing information modules. Will add to later versions...

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


#directory for source files
setwd("C:/Users/u3579238/Work/Software/danrosauer-scripts/sPhy")

##Necessary source files (pls load all)

##code for the wrapper functions that run the landscape simulation and the tree generation
source('C:/Users/u3579238/Dropbox/sDiv S-Phy/groups/simulation/code/simulation/Working_Sims/script_to_run_scape5_simplified.R', chdir = TRUE)

##back end code for the scape script - generates landscapes using trees
source('C:/Users/u3579238/Dropbox/sDiv S-Phy/groups/simulation/code/simulation/Working_Sims/gradient_cart_simulation_v3.1.R', chdir = TRUE)

##plotting and landscape perturbation functions
source('C:/Users/u3579238/Dropbox/sDiv S-Phy/groups/simulation/code/simulation/Working_Sims/landscape_manipulation_plotting.R', chdir = TRUE)

##generate trees with a variation of shapes
source('C:/Users/u3579238/Dropbox/sDiv S-Phy/groups/simulation/code/simulation/Working_Sims/tree_generate_fct.R', chdir = TRUE)

#Meaning of scape types - 1st number is the input # for the runscape function, second # is the figure reference for each landscape type.
	 
	##Environment covary w phylogeny
	#1 = #1a) ##lambda E = 1, lambda R = 0, repulsion
	#weak attraction (i.e. env filtering), strong repulsion(competition), with no signal of range size

	#2 = #1b) ##lambda E = 1, lambda R = 0, attraction
	#strong attraction (i.e. env filtering), weak repulsion (competition), with no signal of range size

	##Environment and range covary w phylogeny
	#3 = #2a)##lambda E = 1, lambda R = 1, repulsion
	
	#4 = #2b)##lambda E = 1, lambda R = 1, attraction
	
	##No environment or range covariance with phylogeny
	#5 = #3 spatial)##lambda E =0, lambda R = 0
	
	#6 = #3 aspatial)##lambda E =0, lambda R = 0
	
	##Range phylogeny covariance with no environmental
	#signal
	#7 = #4 spatial)  lambda E =0, lambda R = 1
	
	#8 = #4 aspatial)##lambda E =0, lambda R = 1,


###############USING SCAPE SETUP################################################
####################Generate trees - creates 10 times as many trees as you ask for, then subsamples in order to increase variation.

numsp= 64 	#specify number of species, must be power of 2
numtrees=1 	#specify number of trees to generate
tree_gen<-tree_generate(numsp,numtrees) 	#generate trees
treeinfo<-tree_gen$treeinfo 	#save resulting tree info
treelist<-tree_gen$treelist		#save resulting trees



###################Generate scapes
#uses default function values, allows you to change tree size, landscape size, specify landscape types to create

scape.size=100 	#specify landscape size (this is one edge of square matrix)
total.cells=(scape.size+1)^2 
maxabund=100 #what is the highest abundance (carrying capacity) possible in a patch?
#subsample.percent=75	
#detect.threshold=5

plotscapes=FALSE #do you want plots of all the landscapes generated? (say no if there are more than 3 or 4)
scape.type=2 	#what scape types? "all" or c(1,2,3,4), etc

##create scapes
Spmatrixlist <- runscape(numsp=numsp,scape.type=scape.type,scape.size=scape.size,treelist,maxabund=maxabund,plotscapes=FALSE) #generate a list of landscapes in the requested scape types using the supplied trees. 

#######################################################
# beyond here code is specific to Marxan preparation  #
#######################################################

library(phylobase)

source("C:/Users/u3579238/Work/Software/danrosauer-scripts/sPhy/Prepare_Marxan_Inputs.r")

# TEMP SET TO A SINGLE SCAPE DURING DEVELOPMENT
scapenum          <- 1
first_dir_num     <- 1
base_dir          <- "C:/Users/u3579238/Work/Software/Marxan/scape_test/"
output_dir        <- "run" 

Spmatrix <- Spmatrixlist[scapenum][[1]]
treenum <- Spmatrix[[3]]
tree   <- phylo4(treelist[[treenum]])

matrix_PA  <- Spmatrix[[2]]

taxa <- dimnames(matrix_PA)[[2]]
taxon_count <- length(taxa)
node_count  <- nrow(tree@edge)
cell_count <- nrow(matrix_PA)
cells <- data.frame(cell = 1:cell_count,area = rep(1,cell_count))

# expand the presence absence matrix to have a column for each internal node
matrix_PA_nodes <- matrix(nrow=cell_count,ncol=node_count)
matrix_PA_nodes[,1:taxon_count] <- matrix_PA

# give names to internal branches (added to species names)
node_labels <- paste("node_", taxon_count+1:nNodes(tree), sep="")
nodeLabels(tree) <- node_labels

nodes <- getNode(tree,1:node_count)
dimnames(matrix_PA_nodes) <- list(NULL,names(nodes))
BranchData  <- data.frame(BranchNames=names(nodes), BranchRange=as.numeric(rep(0,node_count)),stringsAsFactors=FALSE)
BranchData$BranchRange[1:taxon_count] <- apply(matrix_PA_nodes[,1:taxon_count],MARGIN=2,FUN=sum)   

# loop through the internal nodes to add an occurrence column for each to the presence/absence matric
for (m in (taxon_count+1):node_count) {
  
  branch_name  <- names(descendants(tree, m, type="tips"))
  branch_occ    <- apply((matrix_PA[,which(taxa %in% branch_name)]),MARGIN=1,FUN=max)
  matrix_PA_nodes[,m] <- branch_occ
  BranchData$BranchRange[m] <- sum(branch_occ)
  
  node <- as.integer(nodes[m])
  
  # progress output - remove once working
  if (m %% 10 == 0) {
    cat("\n",m,"of",node_count,"done")
    cat("\n",names(nodes)[m],"\t",node_count,"\n")
  }
}

#name the output folder
if (scapenum > 1) {
  directory_number <- scapenum + first_dir_num - 1  # allows for starting with a later number
  this_output_dir <- paste(base_dir,output_dir,"_",directory_number ,sep="")
} else {
  this_output_dir <- paste(base_dir,output_dir,sep="")
}

result <- MarxanInputs(write.dir = this_output_dir,
                       quad_list = cells,
                       tree = tree,
                       branch_data = BranchData,
                       node_occ = matrix_PA_nodes,
                       target =          0.1,
                       spf_multiplier =  10)

gc() # garbage collect from memory
