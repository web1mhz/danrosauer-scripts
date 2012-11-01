
library(SDMTools)
library(doParallel)

geogdist = function(x1,y1,x2,y2) {
  distance = sqrt((x1-x2)^2 + (y1-y2)^2)
  return(distance)
}

####### Parameters  #######
background_samplecount = 10000
s = 1.0  # the standard deviation for the gaussian distance 
cores = 6
###########################


#define directories
work.dir = 'C:/Users/u3579238/Work/Phylofest/U_lithomoda/'; setwd(work.dir)
clipped.bioclim = 'C:/Users/u3579238/Work/Phylofest/U_lithomoda/grids_clipped/'

#define some basic data
bio1.filename <- paste(clipped.bioclim,'bio1_msk.asc',sep="")
samplegrid.filename <- paste(work.dir,'species_all_trimdupes.asc',sep="")

pixels.env = asc2dataframe(bio1.filename)
pixels.sample = asc2dataframe(samplegrid.filename)

background_ind    = sample((1:nrow(pixels.env)),background_samplecount)
pixels.background = pixels.env[background_ind,]
pixels.combined = rbind(pixels.sample,pixels.background)
pixels.combined[,3] = 0
pixels.combined = unique(pixels.combined)

sample.rows = nrow(pixels.sample)
combined.rows = nrow(pixels.combined)

distmatrix = matrix(-1,nrow=sample.rows,ncol=combined.rows)

cl <- makeCluster(cores)
registerDoParallel(cl)

cat("\nAbout to calculate geographic distances on ",cores, " cores\n")
foreach (i = 1:sample.rows) %dopar% {
  for (j in 1:combined.rows) {
    distmatrix[i,j] = geogdist(pixels.sample[i,1],pixels.sample[i,2],pixels.combined[j,1],pixels.combined[j,2])
  }
}

distmatrix = exp ( - (distmatrix ^ 2 ) / ( 2 * s ^ 2 ) )  # gaussian transform of values in distance matrix
pixels.combined[,3] = colSums(distmatrix)                 # 
names(pixels.combined)[3] = "gauss_dist"

cat("\nAbout to merge with full area grid\n")
pixels.env <- merge(pixels.env,pixels.combined,all.x=TRUE)
pixels.env <- pixels.env[c("y","x","gauss_dist")]

cat("\nAbout to write grid to file\n")
dataframe2asc(pixels.env,"bias_grid_s1.asc")

cat("\nFinished\n")
