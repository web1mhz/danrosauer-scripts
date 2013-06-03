# a script to extract environment values for a set of points, from a set of ascii grids 
# Dan Rosauer - October 2012

library(SDMTools)

########## Parameters ##########
#inputs
work.dir          <- 'C:/Users/u3579238/Work/Phylofest/Models/geckoes/sequence_sites'

# expects column 2= latitude, column 3 = longitude
samples.dir       <- 'C:/Users/u3579238/Work/Phylofest/Models/geckoes/sequence_sites/'               
samples.filename  <- 'Phyllurus_lin_loc.csv'
env.dir           <- 'C:/Users/u3579238/GISData/EnvironmentGrids/AusGDMGrids/ForMaxent'

#outputs
output.filename   <- paste(samples.dir,"env_at_",samples.filename,sep="")
################################

setwd(work.dir)

points <- read.csv(paste(samples.dir,samples.filename,sep=""))
#pointsxy <- points[,3:2]
pointsxy <- points[,7:6]

#extract values from each environmental layer in the folder
grids_to_use <- list.files(env.dir,pattern='*.asc',full.name=TRUE)
to_exclude   <- grep("aux.xml",grids_to_use)
to_exclude   <- c(to_exclude, grep("asc.ovr",grids_to_use))
grids_to_use <- grids_to_use[- to_exclude]

for (tfile in grids_to_use) {
  cat("\nabout to do", tfile)
  tasc = read.asc(tfile) #read in the data
  dataname = gsub(env.dir,'',tfile);  dataname = gsub('\\_msk.asc','',dataname)
  dataname = gsub('/','',dataname) #define the column name
  points[dataname] = extract.data(pointsxy,tasc) #append the data
}

cat("\nAbout to write table to file\n")
write.csv(points,output.filename)

cat("\nFinished\n")
