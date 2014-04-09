## this script calculates richness and endemism from modelled suitability surfaces
## it requires all the model grids to have the same extent
rm(list=ls())

library(SDMTools)
library(raster)
library(ape)
library(phylobase)
source("C:/Users/u3579238/Work/Software/danrosauer-scripts/Phylofest_code/phylogenetic endemism.r")

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

max.rows=10000000

#define directories
base.dir   =    'C:/Users/u3579238/Work/AMT/Models/'
#input.dir       <- 'lineage_models/asc_clipped'
input.dir       <- 'lineage_models/asc_clipped_cube_method/'
output.dir      <- base.dir
file_pattern    <- 'lin_model_het'
#file_pattern    <- 'gehyra_'
template_grid   <- 'C:/Users/u3579238/Work/AMT/Models/species_models/maxent/Heteronotia/heteronotia_binoei_17mar14_mean.asc'
#group_lin_file  <- 'group_lineage_list.csv'
group_lin_file  <- 'Heteronotia_binoei_lineage_list.csv'

#tree details  - this works for one genus at a time
#tree.file     = 'trees/TreeLineagesBinoei181113.nex'
tree.file     = 'trees/TopEnd_8Jan14_strictclock_relabel.nex'
outgroup      = 'planiceps'
#outgroup      = 'heteronotia'
preface       = 'lin_model_'

#output_prefix <- "het_"
output_prefix <- "Hb_"
threshold = 0.01  # this is not a species level threshold, but one used for each lineage model

####  end of parameters

setwd(base.dir)
files <- list.files(path = input.dir, pattern = file_pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

setwd(input.dir)
#template.asc = read.asc.gz(files[1])
template.ras = raster(template_grid)
model_rows=nrow(template.ras)
model_cols=ncol(template.ras)

# the original version, excluding NA cells
finite_cells <- which(is.finite(template.ras[]))
pos <- as.data.frame(rowColFromCell(template.ras,finite_cells)) #get all points that have data
rm(finite_cells)

i <- 0

for (tfile in files) {
  checkname=unlist(strsplit(tfile,".",fixed=T))
  if (checkname[length(checkname)]=="asc") {   # only accept filenames ending in .asc
    #tasc = read.asc.gz(tfile)                            #read in the data
    tras = raster(tfile)                                #read in the data    
    newname <- gsub(".asc",'',tfile)
    newname <- tolower(gsub(preface,"",newname))
    pos[newname] <- tras[cbind(pos$row,pos$col)]           #append the data
    pos[(which(pos[newname]< threshold)),newname] <- 0    # set values below the threshold to 0
    pos[(which(is.na(pos[newname]))),newname]     <- 0    # set the nulls to 0    
    i <- i+1
    cat("\n",i,newname,"loaded")
  }
}

cat("\nRemoving unoccupied cells\n")
cat("Before:",nrow(pos),"\n")
rowsums <- apply(pos[,3:ncol(pos)],1,sum,na.rm=T)
pos <- pos[which(rowsums>0),]
rm(rowsums)
cat("After:",nrow(pos),"\n")
max.rows <- min(max.rows,nrow(pos))
gc()

setwd(base.dir)
group_lin_list <-read.csv(group_lin_file)
group_lin_list$lineage <- tolower(group_lin_list$lineage)

# read in the tree
tree_suffix <- unlist(strsplit(tree.file,"[.]"))[2]
if (tree_suffix == "nex") {
  tree <- read.nexus(tree.file)
} else {
  tree <- read.tree(tree.file)
}
if (outgroup != '') {
  tree <- drop.tip(tree,tip=which(tree$tip.label==outgroup))
}
tree <- phylo4(tree)
labels(tree) <- tolower(labels(tree))

# ensure that the tree tips match the model names
model.names <- names(pos)
model.names <- model.names[-(1:2)] # names of all columns except the 1st two which are row, col
#model.names <- tolower(gsub(preface,"",model.names))
model.groups <- data.frame(model.groups=vector("character",nTips(tree)),stringsAsFactors=F)

for (i in 1:nTips(tree)) {
  cat(labels(tree)[i],"\n")
  row <- group_lin_list[group_lin_list$lineage==labels(tree)[i],]
  if (nrow(row) > 0) {
    labels(tree)[i] <- tolower(paste(row$ModelGroup,row$lineage,sep="_"))
    model.groups[i,1] <- as.character(row$ModelGroup)
  }
}

tree <- phylo4d(tree,tip.data=model.groups)

tree.names  <- as.character(labels(tree)[1:nTips(tree)])
matched.names <- intersect(model.names,tree.names)
matched.tips  <- which(labels(tree,"tip") %in% matched.names)
cat("\nNot in tree names:",setdiff(model.names,tree.names),"\n")
cat("\nNot in model names:",setdiff(tree.names,model.names),"\n")

# a subtree containing only the tips for which there is a corresponding model
subtree <- subset(tree,tips.include=matched.tips)

write.csv(pos,paste(output_prefix,"sites_x_lineage.csv",sep=''))

gc()
result <- calc_PE_mymodels(subtree,pos[1:max.rows,which(names(pos) %in% matched.names)])
gc()

pos_output <- cbind(pos[1:max.rows,1:2],result)
pos_output <- pos_output[,-3] # omit the site column

# flip the y values
#pos_output$row <- max(pos_output$row) + 1 - pos_output$row

# add residual columns
PE_WE_mod <- lm(pos_output$PE~pos_output$WE)
pos_output$PE_WE_resid <- PE_WE_mod$residuals
PE_WE_loglog_mod <- lm(log(pos_output$PE)~log(pos_output$WE),subset=which(!is.infinite(log(pos_output$WE))))
pos_output_log <- cbind(pos_output[which(!is.infinite(log(pos_output$WE))),],PE_WE_loglog_mod$residuals)

dataframe2asc(pos_output_log[,c(1,2,8)],paste(output_prefix,"PE_WE_loglog_resid.asc",sep=""),output.dir)

filenames <- c(paste(output_prefix,"PE.asc",sep=""),paste(output_prefix,"PD.asc",sep=""),paste(output_prefix,"WE.asc",sep=""),paste(output_prefix,"SR.asc",sep=""),paste(output_prefix,"PE_WE_resid.asc",sep=""))
dataframe2asc(pos_output[,1:7],filenames,output.dir)

write.csv(pos_output,paste(output_prefix,"scores.csv",sep=""),row.names=FALSE)

# make some output images
PE.ras <- raster(filenames[1])
windows(9,9)
plot(PE.ras,main="PE",col=rainbow(25,start=0.1,end=1))

PD.ras <- raster(filenames[2])
windows(9,9)
plot(PD.ras,main="PD",col=rainbow(25,start=0.1,end=1),ylim=c(180,1400),xlim=c(700,2800))

WE.ras <- raster(filenames[3])
windows(9,9)
plot(WE.ras,main="WE",col=rainbow(25,start=0.1,end=1),ylim=c(180,1400),xlim=c(700,2800))

SR.ras <- raster(filenames[4])
windows(9,9)
plot(SR.ras,main="SR",col=rainbow(25,start=0.1,end=1),ylim=c(180,1400),xlim=c(700,2800))

PE_WE_resid.ras <- raster(filenames[5])
windows(9,9)
plot(PE_WE_resid.ras,main="Residual of PE ~ WE",col=rainbow(25,start=1/6,end=1),ylim=c(240,1400),xlim=c(750,2800))

PE_WE_loglog_resid.ras <- raster(filenames[6])
windows(9,9)
plot(PE_WE_loglog_resid.ras,main="Residual of logPE ~ logWE",col=rainbow(25,start=1/6,end=0.9),ylim=c(240,1400),xlim=c(750,2800))


windows(16,10)
par(mfrow=c(2,2),mar=c(3,4,3,2))
plot(PE.ras,main="PE",col=rainbow(25,start=0.1,end=1),ylim=c(240,1400),xlim=c(750,2800))
plot(PD.ras,main="PD",col=rainbow(25,start=0.1,end=1),ylim=c(240,1400),xlim=c(750,2800))
plot(WE.ras,main="WE",col=rainbow(25,start=0.1,end=1),ylim=c(240,1400),xlim=c(750,2800))
plot(SR.ras,main="SR",col=rainbow(25,start=0.1,end=1),ylim=c(240,1400),xlim=c(750,2800))

# now for some residuals
windows(10,10)
plot(pos_output$WE,pos_output$PE,cex=0.7,xlab="WE",ylab="PE")
PE_WE_mod <- lm(pos_output$PE~pos_output$WE)
abline(PE_WE_mod,col='red')
pos_output$PE_WE_resid <- PE_WE_mod$residuals

windows(10,10)
plot(log(pos_output$WE),log(pos_output$PE),cex=0.7,xlab="WE",ylab="PE")

# PE.tasc <- asc.from.raster(PE.ras)
# image.grid(PE.tasc,filenames[1],zlim=c(0,1),cols=rainbow(10)) #plot the image after logging the actual data