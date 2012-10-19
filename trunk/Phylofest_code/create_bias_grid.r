# a script to create a bias grid for use with Maxent
# Dan Rosauer - October 2012

library(SDMTools)

# define the geographic distance function
geogdist = function(x1,y1,x2,y2) {
  distance = sqrt((x1-x2)^2 + (y1-y2)^2)
  return(distance)
}

########## Parameters ##########
#inputs
work.dir          <- 'C:/Users/u3579238/Work/Phylofest/Diporiphora/'

# expects column 2= latitude, column 3 = longitude
samples.filename  <- 'Diporiphora_all_maxent.csv'                     
env.filename      <- 'C:/Users/u3579238/Work/Phylofest/Diporiphora/Diporiphora_grids_clipped/bio1_msk.asc'

#outputs
output.filename   <- "bias_grid.asc"

#settings
background_samplecount <- 20000
#s <- 100  # the standard deviation for the gaussian distance in kilometres
s <- 2  # the standard deviation for the gaussian distance in degrees
################################

setwd(work.dir)

grid.env      <- read.asc(env.filename)
rowcount.env  <- attributes(grid.env)$dim[2]
colcount.env  <- attributes(grid.env)$dim[1]

# place non-zero values in the corners to ensure output has same extent
corner_ll <- c(1,1)
corner_ul <- c(1,rowcount.env-1)
corner_lr <- c(colcount.env-1,1)
corner_ur <- c(colcount.env-1,rowcount.env-1)
corners   <- rbind(corner_ll,corner_ul,corner_ur,corner_lr)
grid.env[corners] <- 1

# then write the file to disk and read in with corners.  Inefficient, but works!
write.asc(grid.env,"envtemp.asc")

# load pixel coords and values from the environment layer (which defines the model extent)
pixels.env <- asc2dataframe("envtemp.asc")
system("del envtemp.asc")
system("rm envtemp.asc")

#read in sample sites and make a raster
sample_sites <- read.csv(samples.filename)
sample_sites <- sample_sites[c(3,2)]
sample_sites[,3] <- 1

grid.env.att <- attributes(grid.env)
env.xll <- grid.env.att$xll
env.yll <- grid.env.att$yll
env.cellsize <- grid.env.att$cellsize

mat <- matrix(NA,nrow=colcount.env,ncol=rowcount.env)  #strange.  I have reveresed rows and columns in the matrix to get the correct grid shape
grid.sample <- as.asc(mat,xll=env.xll,yll=env.yll,cellsize=env.cellsize)
grid.sample <- put.data(sample_sites,grid.sample)

samplegrid.filename <- paste(work.dir,sub(".csv",".asc",samples.filename),sep="")
write.asc(grid.sample,samplegrid.filename)

pixels.sample <- asc2dataframe(samplegrid.filename)
pixels.sample <- pixels.sample[pixels.sample[3]>0,]
delete_cmd <- paste("del",samplegrid.filename)
system(delete_cmd)

background_ind    <- sample((1:nrow(pixels.env)),background_samplecount)
pixels.background <- pixels.env[background_ind,]
pixels.combined <- rbind(pixels.sample,pixels.background)
pixels.combined[,3] <- 0
pixels.combined <- unique(pixels.combined)

sample.rows <- nrow(pixels.sample)
combined.rows <- nrow(pixels.combined)

distmatrix <- matrix(-1,nrow=sample.rows,ncol=combined.rows)

cat("\nAbout to calculate geographic distances\n")
for (i in 1:sample.rows) {
  for (j in 1:combined.rows) {

    #calculate distance between points in km, allowing for curvature of the earth using Vincenty's formula
    
#     distmatrix[i,j] <- (distance(pixels.sample[i,1],pixels.sample[i,2],pixels.combined[j,1],pixels.combined[j,2])$distance) / 1000

    #calculate distance between points in degrees
    distmatrix[i,j] = geogdist(pixels.sample[i,1],pixels.sample[i,2],pixels.combined[j,1],pixels.combined[j,2])
  }
  if (i%%5==0) {
    cat("\nSample",i,"of",sample.rows)
  }
}

distmatrix <- exp ( - (distmatrix ^ 2 ) / ( 2 * s ^ 2 ) )  # gaussian transform of values in distance matrix
pixels.combined[,3] <- colSums(distmatrix)                 # 
names(pixels.combined)[3] <- "gauss_dist"

cat("\nAbout to merge with full area grid\n")
pixels.env <- merge(pixels.env,pixels.combined,all.x=TRUE)
pixels.env <- pixels.env[c("y","x","gauss_dist")]

cat("\nAbout to write grid to file\n")
dataframe2asc(pixels.env,output.filename)

cat("\nFinished\n")
