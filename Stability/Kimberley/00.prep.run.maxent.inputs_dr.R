#drafted by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
#GNU General Public License .. feel free to use / distribute ... no warranties

################################################################################
rm (list=ls())
library(SDMTools)

################################################################################
################################################################################

#define directories
work.dir            <- 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/'; setwd(work.dir)
current.bioclim     <- 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000'
maxent.jar          <- 'C:/Users/u3579238/Work/Refugia/Stability/maxent.jar'
veg_grid            <- 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/IBRA_NarrowAMT2.asc'
mask_layer_name     <- veg_grid
output_folder_name  <- 'maxent.output_narrow_AMT2'

#define some basic data
mvg.asc  <- read.asc(veg_grid)                 # read in the vegetation grid
area.asc <- mvg.asc; area.asc[which(is.finite(area.asc) & area.asc!=1)] = 0  #set all other areas != 1 to 0

maxent_threads      <- 9

ymin <- -24

################################################################################
#get a subset of the data for occur & background
pos = as.data.frame(which(is.finite(mvg.asc),arr.ind=TRUE)) #get all points that have data
pos$mvg = mvg.asc[cbind(pos$row,pos$col)] #append the area data
pos.subset = pos[which(pos$row %% 2 == 0 & pos$col %% 2 == 0),] #get a subset of the data as a regular grid of every second cell (even row / col numbers)
#pos.subset= pos  #this line is an alternative to the previous, which selects every 2nd cell
pos.subset = rbind(pos.subset,pos[which(pos$mvg==1),]); pos.subset = unique(pos.subset) #ensure all cells in the target area are included in dataset

# apply a southern limit
xy <- getXYcoords(area.asc)
ycoords <- which(xy$y >= ymin)
pos.subset <- pos.subset[pos.subset$row %in% ycoords,]

#append the current environmental data
for (tfile in list.files(current.bioclim,pattern='\\.asc.gz',full.name=TRUE)) {
	tasc = read.asc.gz(tfile) #read in the data
	dataname = gsub(current.bioclim,'',tfile);
  dataname = gsub('\\.asc.gz','',dataname);
  dataname = gsub('/','',dataname)
	pos.subset[dataname] = tasc[cbind(pos.subset$row,pos.subset$col)] #append the data
}
pos.subset = na.omit(pos.subset) #ensure there is no missing data

#now remove the mask column (if any)
if (exists("mask_layer_name")) {
  pos.subset[mask_layer_name] <- NULL
}

#define the occurrences & background ... then write out the data
occur = data.frame(species='area',pos.subset[which(pos.subset$mvg==1),])
occur$mvg = NULL #define the occurrences
bkgd = data.frame(species='bkgd',pos.subset)
bkgd$mvg=NULL #define the background

write.csv(occur,'occur.csv',row.names=FALSE) #write out the occurrences
write.csv(bkgd,'bkgd.csv',row.names=FALSE) #write out the background

################################################################################
#run maxent
dir.create(output_folder_name)
sys_command = paste('java -mx2048m -jar ',maxent.jar,' -e bkgd.csv -s occur.csv -o', output_folder_name,'nothreshold nowarnings novisible -P jackknife -r -a') # variant with jacknifing for better info on variable importance

# to run in parallel
if (maxent_threads > 1) {
  sys_command <- paste(sys_command, " threads=",maxent_threads,sep="")
}

system(sys_command)
