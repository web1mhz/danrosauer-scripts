rm(list=ls())

library(phylobase)
library(ape)

source("C:/Users/u3579238/Work/Software/danrosauer-scripts/Phylofest_code/calc_PE.r")

#define files
species_xy.filename <- "C:/Users/u3579238/Dropbox/Biodiverse presentations/Biodiverse Prac Oct 2013/ModelledHylids8Sep08.csv"
tree.file           <- "C:/Users/u3579238/Work/Study_Taxa/Herps/Hylidae/hylid_relaxed alltaxa 27apr09.nex"
remap.file          <- "C:/Users/u3579238/Work/Study_Taxa/Herps/Hylidae/Translate_Hylid_Names_aug08_2cols.csv"

####  end of parameters

# load the data
species_xy <- read.csv(species_xy.filename)
species_xy$site <- paste(species_xy$Longitude,species_xy$Latitude,sep=":")
tree            <- phylo4(read.nexus(tree.file))
remap           <- read.csv(remap.file,stringsAsFactors = FALSE)

### take a subset of the spatial data, for a quick demo ###
species_xy <- species_xy[(species_xy$Latitude > -25),]
species_xy <- species_xy[(species_xy$Longitude > 140),]

# remap tree names to match the spatial data
labels <- as.character(tipLabels(tree))
for (i in 1:length(labels)) {
  remap.row <- remap[remap$Name_tree==labels[i],]
  if (nrow(remap.row)==1) {tipLabels(tree)[i] <- remap.row$Name_spatial}
}

rm(i,remap,remap.row) # cleaning up

# ensure that the tree tips match the spatial names
spatial_names <- as.character(unique(species_xy$Species))
labels <- as.character(tipLabels(tree))
on_tree     <- intersect(spatial_names,labels)
not_on_tree <- setdiff(spatial_names,labels)

cat("\nNot in tree:\n",not_on_tree,"\n")
cat("\nNot in spatial:\n",setdiff(labels,spatial_names),"\n")

# a subtree containing only the tips for which there is corresponding site data
subtree <- subset(tree,tips.include=on_tree)

#get unique sites
sites_xy <- unique(species_xy[,2:4])

# create sites_x_tips data frame
sites_x_tips <- data.frame(sites_xy$site,row.names=NULL)
names(sites_x_tips) <- "sites"

# create a sites_x_tips matrix with the same order as the tree
for (i in 1:nTips(subtree)) {
  sites_x_tips[,i+1] <- rep(0,nrow(sites_x_tips))
  this_label <- labels(subtree)[i]
  names(sites_x_tips)[i+1] <- this_label
  this_species_xy <- species_xy[which(species_xy$Species==this_label),4]
  sites_x_tips[which(sites_x_tips$site %in% this_species_xy),i+1] <- 1
}

gc()
result <- calc_PE(subtree,sites_x_tips[,-1],"presence") # call calc_PE, omiting the site name column from the matrix
output <- cbind(sites_xy,result) # add on the lat and long columns
gc()


### optional plot PE result ###

library(classInt)
class_count <- 15

my.class.fr<-classIntervals(log(output$PE),n=class_count,style="equal")   
my.col.fr<-findColours(my.class.fr,rainbow(class_count,start=0.1)) # ramp colors

# Map
windows()
plot(output$Longitude,output$Latitude,col=my.col.fr,pch=20)