# a script to prepare site pair, dissimilarity and environment data the GDM tables in R
# Dan Rosauer - October 2012

library(SDMTools)

# define the geographic distance function
geogdist = function(x1,y1,x2,y2) {
  distance = sqrt((x1-x2)^2 + (y1-y2)^2)
  return(distance)
}

########## Parameters ##########
#inputs
work.dir          <- 'C:/Users/Dan/GDM_R/'
env.dir           <- 'C:/Users/Dan/Work/GISData/AsciiGrids/test_set'
site_pair_file    <- 'C:/Users/Dan/Study_Taxa/Herps/Myobatrachidae/phylo_dist_Myob_chrono_0.4_test.csv'
output_file       <- 'sitepairs_env.csv'
################################

setwd(work.dir)

input.sitepairs   <- read.csv(site_pair_file)
new.sitepairs     <- input.sitepairs[c("x0","y0","x1","y1","dist")]

site0_xy <- new.sitepairs[c("x0","y0")]
site1_xy <- new.sitepairs[c("x1","y1")]

#extract values from each environmental layer in the folder
for (tfile in list.files(env.dir,pattern='*.asc',full.name=TRUE)) {
  cat("\nabout to do", tfile)
  tasc = read.asc(tfile) #read in the data
  dataname = gsub(env.dir,'',tfile);  
  dataname = gsub('\\_msk.asc','',dataname)
  dataname = gsub('/','',dataname) #define the column name
  dataname = gsub(".asc",'',dataname);  
  colname0  = paste("site0.",dataname,sep="")
  colname1  = paste("site1.",dataname,sep="")
  new.sitepairs[colname0] = extract.data(site0_xy,tasc) #append the data
  new.sitepairs[colname1] = extract.data(site1_xy,tasc) #append the data  
}

# reorder columns
columns <- c("dist","x0","y0","x1","y1",sort(names(new.sitepairs)[6:ncol(new.sitepairs)]))
new.sitepairs <- new.sitepairs[columns]

new.sitepairs <- na.omit(new.sitepairs)

cat("\nAbout to write table to file\n")
write.csv(new.sitepairs,output_file)

cat("\nFinished\n")
