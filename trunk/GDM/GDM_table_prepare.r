# a script to prepare site pair, dissimilarity and environment data the GDM tables in R
# Dan Rosauer - October 2012 

library(SDMTools)

# define the geographic distance function
geogdist = function(x1,y1,x2,y2) {
  distance = sqrt((x1-x2)^2 + (y1-y2)^2)
  return(distance)
}

########## Parameters ##########
software          <- 'C:/Users/u3579238/Work/GDM/software/GDM4Tables_1_0/GDM_Table_Funcs64_1_0.R'
work.dir          <- 'C:/Users/u3579238/Work/GDM/GDM_table_test/'
env.dir           <- 'C:/Users/u3579238/GISData/EnvironmentGrids/AusGDMGrids/AsciiGrids/'
site_pair_file    <- 'C:/Users/u3579238/Work/Study_Taxa/Herps/Myobatrachidae/species&phylo_Myob_chrono_0.4.csv'
output_file       <- 'sitepairs_env1.csv'
dist.fieldname    <- 'phylo_sorenson'

grids_to_use      <- c("microgi.asc","mesogi.asc","C4gi.asc","rhu215_x.asc","arid_max.asc","rprecmax.asc","slrain1.asc","slrain2.asc","radni.asc","radnx.asc","mintx.asc","rtxmin.asc","distpermwat_limit1.asc","gravity.asc","geollrngeage.asc","nmnln0.asc","ridgetopflat.asc","slope.asc","mvs31_ccrf.asc")
# comment the precding line and get all grids in the directory if grids_to_use doesn't exist
################################

setwd(work.dir)
source(software)

input.sitepairs   <- read.csv(site_pair_file)
names(input.sitepairs)[which(names(input.sitepairs)==dist.fieldname)] <- "observed"
new.sitepairs     <- input.sitepairs[c("observed","x0","y0","x1","y1")]
new.sitepairs["weight"] <- rep.int(1,nrow(new.sitepairs))

site0_xy <- new.sitepairs[c("x0","y0")]
site1_xy <- new.sitepairs[c("x1","y1")]

# this gets all grids in the floder, if none specified
if (!exists("grids_to_use")) {
  grids_to_use <- list.files(env.dir,pattern='*.asc',full.name=TRUE)
}

#extract values from each environmental layer in the folder
for (tfile in grids_to_use) {
  cat("\nabout to do", tfile)
  file_path = paste(env.dir,tfile,sep="")
  tasc = read.asc(file_path) #read in the data
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
columns <- c("observed","weight","x0","y0","x1","y1",sort(names(new.sitepairs)[7:ncol(new.sitepairs)]))
new.sitepairs <- new.sitepairs[columns]

new.sitepairs <- na.omit(new.sitepairs)

below.zero    <- which(new.sitepairs$observed<0)
cat("\n",length(below.zero),"records removed due to a dissimilarity < 0\n")
new.sitepairs <- new.sitepairs[which(new.sitepairs$observed>=0),]


cat("\nAbout to write table to file\n")
write.csv(new.sitepairs,output_file)

rm(input.sitepairs, site0_xy, site1_xy)

cat("\nFinished at",date(),"\n")
