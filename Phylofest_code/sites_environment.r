# a script to extract environment values for a set of points, from a set of ascii grids 
# Dan Rosauer - October 2012

library(SDMTools)

# define the geographic distance function
geogdist = function(x1,y1,x2,y2) {
  distance = sqrt((x1-x2)^2 + (y1-y2)^2)
  return(distance)
}

########## Parameters ##########
#inputs
work.dir          <- 'C:/Users/u3579238/GISData/Helping/Sally'

# expects column 2= latitude, column 3 = longitude
samples.filename  <- 'MaxEntmodel_Wyulda.csv'                     
env.dir      <- 'C:/Users/u3579238/GISData/Helping/Sally/Env_grids_clipped'

#outputs
output.filename   <- paste("env_at_",samples.filename,sep="")
################################

setwd(work.dir)

points <- read.csv(samples.filename)
pointsxy <- points[,3:2]

#extract values from each environmental layer in the folder
for (tfile in list.files(env.dir,pattern='*.asc',full.name=TRUE)) {
  cat("\nabout to do", tfile)
  tasc = read.asc(tfile) #read in the data
  dataname = gsub(env.dir,'',tfile);  dataname = gsub('\\_msk.asc','',dataname)
  dataname = gsub('/','',dataname) #define the column name
  points[dataname] = extract.data(pointsxy,tasc) #append the data
}

cat("\nAbout to write table to file\n")
write.csv(points,output.filename)

cat("\nFinished\n")
