
library(SDMTools)

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

generate_coords = function(rowcount,colcount) {
  rows= c()
  cols= c()
  
  #generate rows
  for (r in 1:rowcount) {
    rr=rep(r,colcount)
    rows = c(rows,rr)
    
    cc=1:colcount
    cols=c(cols,cc)
  }
  
  rc = cbind(rows,cols)
  return(rc)
}
################################################################################
################################################################################

#define directories
#base.dir   = 'C:/Users/u3579238/Work/Refugia/DistModels_Reside/Reptiles_cropped/'; setwd(base.dir)
base.dir   = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models_ext/'; setwd(base.dir)
output.dir = 'C:/Users/u3579238/Work/Refugia/Results/'
file_suffix = ".asc"  # frog models end in _1999.asc.gz, reptile models don't

richness_output = "rept_rich_lin" #"frog_rich_sp_thr0.5"
endemism_output = "rept_end_lin"  #"frog_end_sp_thr0.5"

threshold = 0.0  # this is not a species level threshold, but a generic threshold across all the models

files = list.files(path = base.dir, pattern = file_suffix, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

#first_model.asc = read.asc.gz(files[1])
first_model.asc = read.asc(files[1])
model_rows=nrow(first_model.asc)
model_cols=ncol(first_model.asc)

# the original version, excluding NA cells
pos = as.data.frame(which(is.finite(first_model.asc),arr.ind=TRUE)) #get all points that have data

# trying a list of coords based on the grid size
#pos = data.frame(generate_coords(model_rows,model_cols))

for (tfile in files) {
  checkname=unlist(strsplit(tfile,".",fixed=T))
  if (checkname[length(checkname)]=="asc") {   # only accept filenames ending in .asc
    #tasc = read.asc.gz(tfile)                            #read in the data
    tasc = read.asc(tfile)                                #read in the data    
    dataname=gsub(file_suffix,'',tfile)
    if (dataname != "Cophixalus_peninsularis") {              #skipping a dodgy model

      pos[dataname] = tasc[cbind(pos$row,pos$col)]            #append the data
      pos[(which(pos[dataname]< threshold)),dataname] <- 0    # set values below the threshold to 0
      pos[(which(is.na(pos[dataname]))),dataname]     <- 0    # set the nulls to 0    
      cat("\n",dataname," loaded")
    }
  }
}

cat("\nData loaded, now performing calculations")

# add the suitability values across all models
cols = ncol(pos)
pos$total = apply(pos[3:cols],1,'sum')      # analogous to richness

# the endemism calculation
species_sums = apply(pos[3:cols],2,'sum')
pos_end = pos[1:cols]
for (species in names(pos[3:cols])) {
  if (species_sums[species] > 0) {    # avoid division by zero if a species column sums to 0
    pos_end[species] <- (pos_end[species] / species_sums[species])
  }
}

pos$we = apply(pos_end[3:cols],1,'sum')
rm(pos_end)

pos_output <- pos[,c(2,1,cols+1,cols+2)]

dataframe2asc(pos_output,c(richness_output,endemism_output),output.dir)

# write.asc(rich.asc, paste(output.dir,'sp_richness.asc',sep='') #write out the ascii grid file
# image.grid(rich.asc,paste(output.dir,'sp_richness.png',sep=''),zlim=c(0,1),cols=model_cols) #plot the image after logging the actual data