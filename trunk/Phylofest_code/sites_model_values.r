# a script to extract environment values for a set of points, from a set of ascii grids 
# Dan Rosauer - October 2012

library(SDMTools)

########## Parameters ##########
#inputs
work.dir          <- 'C:/Users/u3579238/Work/Phylofest/Models/skinks/L_delicata_Tingley/'
samples.dir       <- 'C:/Users/u3579238/Work/Phylofest/Models/skinks/L_delicata_Tingley/'              
samples.filename  <- 'Ldelicata_ALA.csv'
lat_col           <- 1
lon_col           <- 2
env.dir           <- 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models/'
env.pattern       <- 'lin_model_lampropholis_delicata_tingley_dr(.)+.asc'  #regex
minimum_value     <- 0.0005

#outputs
name = substr(samples.filename,1, nchar(samples.filename)-4)
output.filename   <- paste(samples.dir,"clades_at_",name,"_dr.csv",sep="")
################################

setwd(work.dir)

points <- read.csv(paste(samples.dir,samples.filename,sep=""))
pointsxy <- points[,c(lon_col,lat_col)]

#extract values from each environmental layer in the folder
grids_to_use <- list.files(env.dir,pattern=env.pattern,full.names=TRUE)
to_exclude   <- grep("aux.xml",grids_to_use)
to_exclude   <- c(to_exclude, grep("asc.ovr",grids_to_use))
#to_exclude   <- c(to_exclude, grep("tingley_dr",grids_to_use))
grids_to_use <- grids_to_use[- to_exclude]

for (tfile in grids_to_use) {
  cat("\nabout to do", tfile)
  tasc = read.asc(tfile) #read in the data
  dataname = gsub(env.dir,'',tfile);  dataname = gsub('\\_msk.asc','',dataname)
  dataname = gsub('/','',dataname) #define the column name
  points[dataname] = round(extract.data(pointsxy,tasc),4) #append the data
  points[dataname][points[dataname]<minimum_value] <- 0  # set values below minimum_value to 0
}

cat("\nAbout to write table to file\n")
write.csv(points,output.filename)

cat("\nFinished\n")
