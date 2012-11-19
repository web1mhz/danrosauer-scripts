# a script to create a bias grid for use with Maxent
# Dan Rosauer - November 2012

library(SDMTools)
library(doParallel)

# define the geographic distance function
geogdist = function(x1,y1,x2,y2) {
  distance = sqrt((x1-x2)^2 + (y1-y2)^2)
  return(distance)
}

########## Parameters ##########
#inputs
taxon_name        <- "Hypsilurus"
higher_taxon = "dragons"
base.dir = paste("C:/Users/u3579238/work/Phylofest/Models/",higher_taxon,"/",sep="")
samples.filename = paste("species_sites/",taxon_name,"_maxent.csv",sep="")
lat_column=2
long_column=3

#settings
background_samplecount <- 50000
#s <- 100  # the standard deviation for the gaussian distance in kilometres
s <- 3  # the standard deviation for the gaussian distance in degrees

env.filename      <- paste("species_models/clipped_grids/",taxon_name,"/bio01_msk.asc",sep="")
#output.filename   <- paste("species_models/bias_grids/",taxon_name,"/bias_grid_test.asc",sep="")
output.filename   <- paste("/",taxon_name,"_bias_grid_",s,"deg.asc",sep="")
cores = 2
################################

work.dir       <- paste(base.dir,"species_models/bias_files",sep="")
setwd(work.dir)
command = paste("cd",work.dir)
system(command)

env.filename  <- paste(base.dir,env.filename,sep="")
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
samples.file <- paste(base.dir,samples.filename,sep="")
sample_sites <- read.csv(samples.file)
#sample_sites <- sample_sites[sample_sites$Use == 1,]

# get the grid attributes
grid.env.att <- attributes(grid.env)
env.cellsize <- grid.env.att$cellsize
# note, these corner coordinates are for the centre of the corner cells, so 1/2 a cell in from the actual grid edges
env.xll <- grid.env.att$xll  
env.yll <- grid.env.att$yll
env.xur <- env.xll + corner_ur[1] * env.cellsize
env.yur <- env.yll + corner_ur[2] * env.cellsize

#remove sample sites outside the boundaries of grid.env
sample_sites <- sample_sites[sample_sites[long_column]>=(env.xll - (env.cellsize * 0.5)),]
sample_sites <- sample_sites[sample_sites[long_column]<=(env.xur + (env.cellsize * 0.5)),]
sample_sites <- sample_sites[sample_sites[lat_column]>=(env.yll - (env.cellsize * 0.5)),]
sample_sites <- sample_sites[sample_sites[lat_column]<=(env.yur + (env.cellsize * 0.5)),]
sample_sites <- sample_sites[,c(long_column,lat_column)]
sample_sites[,3] <- 1

mat <- matrix(NA,nrow=colcount.env,ncol=rowcount.env)  #  strange.  I have reveresed rows and columns in the matrix to get the correct grid shape
grid.sample <- as.asc(mat,xll=env.xll,yll=env.yll,cellsize=env.cellsize)
grid.sample <- put.data(sample_sites,grid.sample)

samplegrid.filename <- paste(taxon_name,"_sites_grid_temp.asc",sep="")
write.asc(grid.sample,samplegrid.filename)

pixels.sample <- asc2dataframe(samplegrid.filename)
pixels.sample <- pixels.sample[pixels.sample[3]>0,]
delete_cmd <- shQuote(paste("del",samplegrid.filename))
system(delete_cmd)

background_ind    <- sample((1:nrow(pixels.env)),background_samplecount)
pixels.background <- pixels.env[background_ind,]
pixels.combined <- rbind(pixels.sample,pixels.background)
pixels.combined[,3] <- 0
pixels.combined <- unique(pixels.combined)

sample.rows <- nrow(pixels.sample)
combined.rows <- nrow(pixels.combined)

#distmatrix <- matrix(-1,nrow=sample.rows,ncol=combined.rows)
distrow    <- vector("numeric",combined.rows)

cl <- makeCluster(cores)
registerDoParallel(cl)

cat("\nAbout to calculate geographic distances\n")
distmatrix <- foreach (i = 1:sample.rows, .combine=rbind, .verbose=TRUE) %dopar% {
  for (j in 1:combined.rows) {

    #calculate distance between points in km, allowing for curvature of the earth using Vincenty's formula
    
    # distmatrix[i,j] <- (distance(pixels.sample[i,1],pixels.sample[i,2],pixels.combined[j,1],pixels.combined[j,2])$distance) / 1000

    #calculate distance between points in degrees
    distrow[j] = geogdist(pixels.sample[i,1],pixels.sample[i,2],pixels.combined[j,1],pixels.combined[j,2])
    }
    if (i%%5==0) {
      cat("\nSample",i,"of",sample.rows)
    }
  distrow
}

stopCluster(cl)

distmatrix <- exp ( - (distmatrix ^ 2 ) / ( 2 * s ^ 2 ) )  # gaussian transform of values in distance matrix
pixels.combined[,3] <- colSums(distmatrix)                 # 
names(pixels.combined)[3] <- "gauss_dist"

cat("\nAbout to merge with full area grid\n")
pixels.env <- merge(pixels.env,pixels.combined,all.x=TRUE)
pixels.env <- pixels.env[c("y","x","gauss_dist")]

cat("\nAbout to write grid to file\n")
#output.filename = paste(base.dir,output.filename,sep="")
dataframe2asc(pixels.env,output.filename)

cat("\nFinished.\n",date(),"\n")
