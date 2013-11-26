## this script calculates richness and endemism from modelled suitability surfaces
## it requires all the model grids to have the same extent
rm(list=ls())

library(SDMTools, Raster)
library(phylobase)
source("C:/Users/Dan/Work/Software/Phylofest_code/phylogenetic endemism.r")

################################################################################
#first define some functions

#create images  (create images function from Jeremy Vanderwall stability script)
image.grid = function(tasc,tfile,zlim=NULL,cols=NULL) {
  if (is.null(cols)) { cols = colorRampPalette(c('green','brown','yellow','orange','red'))(101) }
  #plot the image
  png(tfile,height=847/100,width=1080/100,units='cm',res=300,pointsize=5)
  layout(matrix(c(rep(rep(1,8),6),1,2,2,2,2,1,1,1),nr=7,nc=8,byrow=TRUE))
  par(mar=c(0,0,0,0),cex=1)
  if (is.null(zlim)) { zlim = range(as.vector(tasc),na.rm=T) }
  image(tasc,ann=FALSE,axes=FALSE,zlim=zlim,col=cols)
  legend.gradient(cbind(c(113,114.5,114.5,113),c(-44,-44,-38,-38)),limits=round(zlim,2),title='stability',cols=cols,cex=1.2)  
  plot(Oz.shape,add=T,border="black",pbg="transparent",lwd=0.6) #add the subregions		
  dev.off()
}

################################################################################
################################################################################

#TEMP
row.limit=100000
#TEMP

#define directories
base.dir   =    'C:/Users/Dan/Work/AMT/Models/'
input.dir       <- 'lineage_models/asc_clipped'
output.dir      <- base.dir
file_suffix     <- '.asc'
template_grid   <- 'C:/Users/Dan/Work/AMT/Models/AMT_template.asc'
group_lin_list  <- 'group_lineage_list.csv'

#tree details  - this works for one genus at a time
tree.file     = 'trees/Gehyra uniq lineages 141113 final.tre'
outgroup      = 'Heteronotia'
preface       = 'lin_model_'

#richness_output = "reptfrog_rich_lin_25Sep_thresh_01"
endemism_output = "gehyra_PE"  #"frog_end_sp_thr0.5"
PE_output <- "gehyra_PE"
threshold = 0.01  # this is not a species level threshold, but one used for each lineage model

####  end of parameters

setwd(base.dir)
files = list.files(path = input.dir, pattern = file_suffix, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

setwd(input.dir)
#template.asc = read.asc.gz(files[1])
template.asc = read.asc(template_grid)
model_rows=nrow(template.asc)
model_cols=ncol(template.asc)

# the original version, excluding NA cells
pos <- as.data.frame(which(is.finite(template.asc),arr.ind=TRUE)) #get all points that have data

i <- 0

for (tfile in files) {
  checkname=unlist(strsplit(tfile,".",fixed=T))
  if (checkname[length(checkname)]=="asc") {   # only accept filenames ending in .asc
    #tasc = read.asc.gz(tfile)                            #read in the data
    tasc = read.asc(tfile)                                #read in the data    
    dataname=gsub(file_suffix,'',tfile)
    if (dataname != "Cophixalus_peninsularis") {              #skipping a dodgy model
      newname <- tolower(gsub(preface,"",dataname))
      pos[newname] <- tasc[cbind(pos$row,pos$col)]           #append the data
      pos[(which(pos[newname]< threshold)),newname] <- 0    # set values below the threshold to 0
      pos[(which(is.na(pos[newname]))),newname]     <- 0    # set the nulls to 0    
      i <- i+1
      cat("\n",i,newname,"loaded")
    }   
  }
}

rowsums <- apply(pos[,3:ncol(pos)],1,sum,na.rm=T)
pos <- pos[which(rowsums>0),]
rm(rowsums)
gc()

setwd(base.dir)
group_lin_list <-read.csv(group_lin_list)

# read in the tree
tree <- read.tree(tree.file)
tree <- drop.tip(tree,tip=which(tree$tip.label==outgroup))
tree <- phylo4(tree)

# ensure that the tree tips match the model names
model.names <- names(pos[3:(i+2)])  # names of all columns except the 1st two which are row, col
#model.names <- tolower(gsub(preface,"",model.names))
model.groups <- vector("character",1)

for (i in 1:nTips(tree)) {
  row <- group_lin_list[group_lin_list$lineage==labels(tree)[i],]
  labels(tree)[i] <- tolower(paste(row$ModelGroup,row$lineage,sep="_"))
  model.groups[i] <- as.character(row$ModelGroup)
}

tree.names  <- as.character(labels(tree)[1:nTips(tree)])
matched.names <- intersect(model.names,tree.names)
cat("Not in tree names:",setdiff(model.names,tree.names))
cat("Not in model names:",setdiff(tree.names,model.names))

gc()
#result <- calc_PE(tree,pos[1:row.limit,which(names(pos) %in% tree.names)],"probability")
#result <- calc_PE{(tree,pos[,which(names(pos) %in% tree.names)],"probability")
result <- calc_PE_mymodels(tree,pos[1:row.limit,which(names(pos) %in% tree.names)],model.groups)
gc()

# create a data frame of branches with
  # a) branch name - equalling column name in pos for the tips
  # b) branch length
  # c) model group for tip
        # and a list of model groups for internal branches?

# loop through the non-tip branches, and for each
  # add a column to pos
  # identify all the descendent tips
  # sum the values within each tip (for each row = grid cell) 
  # combine the values between model groups using:
        # 1 - sum(1 - model value)
  # this is the branch model value

# scale the value of each column to sum to its branch length (instead of to 1 for species model)

# row sums are analogous to PD
# PE is the sum across the row of cell/column total

# cat("\nData loaded, now performing calculations")
# 
# # add the suitability values across all models
# cols = ncol(pos)
# pos$total = apply(pos[3:cols],1,'sum')      # analogous to richness
# 
# # the endemism calculation
# species_sums = apply(pos[3:cols],2,'sum')
# pos_end = pos[1:cols]
# for (species in names(pos[3:cols])) {
#   if (species_sums[species] > 0) {    # avoid division by zero if a species column sums to 0
#     pos_end[species] <- (pos_end[species] / species_sums[species])
#   }
# }
# 
# pos$we = apply(pos_end[3:cols],1,'sum')
# rm(pos_end)

pos_output <- cbind(pos[1:row.limit,2:1],result$PE)
#pos_output <- pos[,c(2,1,cols+1,cols+2)]

dataframe2asc(pos_output,PE_output,output.dir)

result.ras <- raster(paste(output.dir,PE_output,".asc",sep=""))
windows(9,9)
plot(result.ras)

#write.asc(.asc, paste(output.dir,'sp_richness.asc',sep='') #write out the ascii grid file
image.grid("PE_output.asc",paste(output.dir,'PE_gehyra.png',sep=''),zlim=c(0,1),cols=model_cols) #plot the image after logging the actual data