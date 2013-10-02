# prepare Global Mammal Phylogenetic input files for Marxan
# Dan Rosauer - started September 2013

rm(list=ls())
library(gdata)

read.dir              <- "C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/"
write.dir             <- "C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/input/"
branch_lengths.file   <- "mammal_branch_lengths.csv"
branch_ranges.file    <- "mammal_branch_ranges.csv"
node_by_cell.file     <- "mammal_node_by_cell_bd_export.csv"
target                <- 0.1
spf_multiplier        <- 3  # adjust the importance of meeting all conservation targets (relative to cost)

############################################

fixnames <- function (text_in) {
#   if (grep("___",text_in)==1) {
#     text_in <- paste("b",text_in,sep="")
#   }
  text_out <- gsub('___',"",text_in)
  text_out <- gsub(" ","_",text_out)
  text_out <- paste("x_",text_out,sep="")
  return(text_out)
}

############################################
setwd(read.dir)

branch_lengths <- read.csv(branch_lengths.file)
branch_lengths$Key <- fixnames(branch_lengths$Key) #remove spaces and trailing underscores
branch_lengths <- branch_lengths[order(branch_lengths$Key),]

cons.feature.count <- nrow(branch_lengths)
branch_ranges <- read.csv(branch_ranges.file)
branch_ranges$node <- fixnames(branch_ranges$node)
branch_ranges <- branch_ranges[order(branch_ranges$node),]

branch_ranges$node <- as.character(branch_ranges$node)
branch_ranges$target <- ceiling(branch_ranges$range * target)

spec <- data.frame(cbind(1:cons.feature.count,branch_ranges$target,branch_lengths$Value * spf_multiplier,branch_ranges$node))
names(spec) <- c("species_id","target","spf","name")

node_by_cell <- read.csv(node_by_cell.file)
node_by_cell$Key <- fixnames(node_by_cell$Key)
pu_list <- data.frame(unique(node_by_cell[,c("ELEMENT","Axis_0","Axis_1")]))
pu_list <- pu_list[order(pu_list$ELEMENT),]
pu_list$ELEMENT <- as.character(pu_list$ELEMENT)
pu_list <- data.frame(cbind(rep(1:length(pu_list[,1])),pu_list))
names(pu_list) <- c("pu_id","pu_name","Axis_0","Axis_1")
pu.count <- length(pu_list[,1])
pu <- data.frame(pu_list$pu_id)

node_by_cell_SpID <- merge(spec, node_by_cell, by.x="name", by.y="Key")
puvspr2 <- data.frame(node_by_cell_SpID$species_id,node_by_cell_SpID$ELEMENT,rep(1,nrow(node_by_cell)))
names(puvspr2) <- c("species","pu_name","amount")  # this data frame has the species number, but the PU name
rm(node_by_cell_SpID)

puvspr2 <- merge(pu_list, puvspr2, by="pu_name")
puvspr2 <- puvspr2[,c(5,2,6)]
names(puvspr2) <- c("species","pu","amount")

# ensure correct column names for maxent format
names(spec)[1] <- "id"
names(pu) <- "id"

setwd(write.dir)
write.fwf(spec,"spec.dat")
write.fwf(pu,"pu.dat")
write.fwf(puvspr2,"puvspr2.dat")
write.csv(pu_list,"pu_id_lookup.csv")