## this script calculates richness and endemism from modelled suitability surfaces
## it requires all the model grids to have the same extent

library(phylobase)

parent.prob <- function(probabilities) {
  # probabilities is a vector of values between 0 and 1
  # add code to check values are of correct type!
  parent.prob <- 1 - prod(1-probabilities)
  return(parent.prob)
}

scale.to <- function(vec,vec.sum) {
  #mat is a vector
  #this function rescales each the vector values to sum to 'sum'
  vec.tot <- sum(vec,na.rm=TRUE)
  if (vec.tot > 0) {
    vec.out <- vec.sum*vec/vec.tot
  } else {
    vec.out <- rep(0,times=length(vec))  #should columns for tips / branches with no occurrences be removed?
  }
  return(vec.out)
}

calc_PE <- function(tree, sites_x_tips,presence=c("presence","abundance","probability")) {
  # add code to check that the values are correct for the presence type:
    # 0 or 1 for presence - this calculates PE (Rosauer et al 2009)
    # from 0 to 1 for probability - this calculates model weighted PE (Rosauer, in prep)
    # any value for abundance - this calculation is equivalent to BED (Cadotte & Davies 2010)
  
  # change to a phylobase phylo4 object
  if (class(tree) == "phylo") {tree <- phylo4(tree)}
  
  sites_x_branches <- data.frame(cbind(rep(0,nrow(sites_x_tips))))

  for (i in 1:nTips(tree)) {
    sites_x_branches[,i] <- sites_x_tips[,which(labels(tree)[i]==names(sites_x_tips))]
    names( sites_x_branches)[i] <- labels(tree)[i]
    cat(i,dim(sites_x_branches),"\n")
  }
  rm(sites_x_tips); gc()
  branch_labels <- labels(tree)
  branchcount <- length(labels(tree))

  # add names and occupancy columns for internal branches
  for (i in (nTips(tree)+1):branchcount) {
    branch_labels[i] <- paste("b",i,sep="")
    desc <- as.integer(descendants(tree,i, type="tips"))
    if (presence=="abundance") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=sum))
    } else if (presence=="presence") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=max))
    } else if (presence=="probability") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=parent.prob))
    }
    sites_x_branches[,i] <- branch_col
    names(sites_x_branches[i]) <- branch_labels[i]
    cat(i,branch_labels[i],length(desc),"\n")
    gc(verbose=F)
  }
  
  #scale columns (branches) to sum to 1
  sites_x_branches <- apply(sites_x_branches,MARGIN=2,FUN=scale.to,1)
  
  branch.lengths <- as.numeric(edgeLength(tree,1:branchcount))
  sites_x_branches <- sites_x_branches[,1:branchcount] * branch.lengths
  PE.vec <- apply(sites_x_branches,MARGIN=1,FUN=sum,na.rm=T)
  
  PE <- data.frame(cbind(1:nrow(sites_x_branches),PE.vec))
  names(PE) <- c("site","PE")
  return(PE)
}
