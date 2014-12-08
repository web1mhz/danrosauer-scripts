library(phylobase)
  
source("C:/Users/u3579238/Work/Software/danrosauer-scripts/Marxan/Prepare_Marxan_Inputs.r")


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

#######################################
# test for first tree, first matrix #
treenum <- 1
matnum  <- 1
run_count <- 1
j <- 1
first_dir_num <- 1
base_dir <- "C:/Users/u3579238/Work/Software/Marxan/scape_test/"
output_dir <- "run" 
#######################################

tree   <- phylo4(treelist[[treenum]])
matrix_obj <- Spmatrixlist[matnum]
matrix_PA  <- matrix_obj[[1]][[2]]

taxa <- dimnames(matrix_PA)[[2]]
taxon_count <- length(taxa)
node_count  <- nrow(tree@edge)

# the quad list includes the name and area of each area.
#   it may not be needed for use with scape
cell_count <- nrow(matrix_PA)
quads <- data.frame(cell = 1:cell_count,area = rep(1,cell_count))

# expand the presence absence matrix to have a column for each internal node
matrix_PA_nodes <- matrix(nrow=cell_count,ncol=node_count)
matrix_PA_nodes[,1:taxon_count] <- matrix_PA

# give names to internal branches (added to species names)
node_labels <- paste("node_", 1:nNodes(tree), sep="")
nodeLabels(tree) <- node_labels

nodes <- getNode(tree,1:node_count)
dimnames(matrix_PA_nodes) <- list(NULL,names(nodes))
BranchData  <- data.frame(BranchNames=names(nodes), BranchRange=as.numeric(rep(0,node_count)),stringsAsFactors=FALSE)

BranchData$BranchRange[1:taxon_count] <- apply(matrix_PA_nodes,MARGIN=2,FUN=sum)   

branch_occ    <- apply((matrix_PA[,which(taxa %in% branch_name)]),MARGIN=1,FUN=max)
matrix_PA_nodes[,m] <- branch_occ
BranchData$BranchRange[m] <- sum(branch_occ)

tree_table <- as.data.frame(print(tree))
is.tip <- tree_table$node.type == "tip"

for (m in (taxon_count+1):node_count) {

  branch_name  <- names(descendants(tree, m, type="tips"))
  branch_occ    <- apply((matrix_PA[,which(taxa %in% branch_name)]),MARGIN=1,FUN=max)
  matrix_PA_nodes[,m] <- branch_occ
  BranchData$BranchRange[m] <- sum(branch_occ)
  
  node <- as.integer(nodes[m])
  
  # progress output - remove once working
  if (m %% 25 == 0) {
    cat("\n",m,"of",node_count,"done")
    cat("\n",tree_table$label[m],"\t",node_count,"\n")
  }
}

#name the output folder
if (run_count > 1) {
  directory_number <- j + first_dir_num - 1  # allows for starting with a later number
  this_output_dir <- paste(base_dir,output_dir,"_",directory_number ,sep="")
} else {
  this_output_dir <- paste(base_dir,output_dir,sep="")
}

if (output_Type == "Marxan") {
  result <- MarxanInputs(write.dir = this_output_dir,
                         quad_list = quads,
                         tree = tree,
                         branch_data = BranchData,
                         node_occ = matrix_PA_nodes,
                         target =          'func',
                         spf_multiplier =  1)
}

gc() # garbage collect from memory

