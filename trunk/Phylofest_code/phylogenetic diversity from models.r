## this script calculates richness and endemism from modelled suitability surfaces
## it requires all the model grids to have the same extent
rm(list=ls())

library(SDMTools)
library(raster)
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
#row.limit=20000
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

output_prefix <- "gehyra_"
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

cat("\nRemoving unoccupied cells\n")
cat("Before:",nrow(pos),"\n")
rowsums <- apply(pos[,3:ncol(pos)],1,sum,na.rm=T)
pos <- pos[which(rowsums>0),]
rm(rowsums)
cat("After:",nrow(pos),"\n")
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
#result <- calc_PE_mymodels(tree,pos[1:row.limit,which(names(pos) %in% tree.names)],model.groups)
result <- calc_PE_mymodels(tree,pos[,which(names(pos) %in% tree.names)],model.groups)
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

pos_output <- cbind(pos[,2:1],result)
pos_output <- pos_output[,-3] # omit the site column

# add residual columns
PE_WE_mod <- lm(pos_output$PE~pos_output$WE)
pos_output$PE_WE_resid <- PE_WE_mod$residuals
PE_WE_loglog_mod <- lm(log(pos_output$PE)~log(pos_output$WE),subset=which(!is.infinite(log(pos_output$WE))))
pos_output_log <- cbind(pos_output[which(!is.infinite(log(pos_output$WE))),],PE_WE_loglog_mod$residuals)
dataframe2asc(pos_output_log[,c(1,2,8)],paste(output_prefix,"PE_WE_loglog_resid.asc",sep=""),output.dir)



filenames <- c(paste(output_prefix,"PE.asc",sep=""),paste(output_prefix,"PD.asc",sep=""),paste(output_prefix,"WE.asc",sep=""),paste(output_prefix,"SR.asc",sep=""),paste(output_prefix,"PE_WE_resid.asc",sep=""),paste(output_prefix,"PE_WE_loglog_resid.asc",sep=""))
dataframe2asc(pos_output,filenames,output.dir)

write.csv(pos_output,"gehyra_scores.csv")

# make some output images
PE.ras <- raster(filenames[1])
windows(9,9)
plot(PE.ras,main="PE",col=rainbow(25,start=0.1,end=1),ylim=c(180,1400),xlim=c(700,2800))

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


PE.tasc <- asc.from.raster(PE.ras)
image.grid(PE.tasc,filenames[1],zlim=c(0,1),cols=rainbow(10)) #plot the image after logging the actual data