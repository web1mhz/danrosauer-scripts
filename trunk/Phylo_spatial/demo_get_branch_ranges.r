
rm(list=ls())

library(ape)
source("get_branch_ranges.r")

# parameters
n <- 25   # number of species
m <- 100  # number of grid cells

tree <- rcoal(n)  # generate a random (coalescent) tree with n tips

species_names  <- paste("sp",1:n,sep="")
tree$tip.label <- species_names #name the tips to match the species

# populate the cells with species
cells       <- 1:m
cells_x_sp   <- data.frame(cell_id = "", sp_id = "")[- 1,]  # create an empty data frame with 0 rows

for (i in 1:n) {  # loop through the species, assigning them to cells
  range       <- as.integer(runif(1,0,m))
  sp_cells    <- sample(x=cells,size=range,replace=FALSE)
  new_rows    <- cbind(sp_cells,rep(i,range))
  cells_x_sp  <- rbind(cells_x_sp,new_rows)
}

names(cells_x_sp) <- c("cell_id","sp_id")

