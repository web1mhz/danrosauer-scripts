#drafted by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
#GNU General Public License .. feel free to use / distribute ... no warranties

################################################################################
library(SDMTools)

################################################################################
################################################################################

#define directories
work.dir =        'C:/Users/u3579238/GISData/Helping/Sally/Wyulda_stability/'; setwd(work.dir)
current.bioclim = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000'
maxent.jar =      'C:/Users/u3579238/Work/Refugia/Stability/maxent.jar'
occ_grid = 'C:/Users/u3579238/GISData/Helping/Sally/Wyulda_stability/wyulda_sites_e.asc'
species_name="Wyulda"
threads=9

output_folder_name = 'maxent.output_2_5deg_topo'
mask_layer_name='wyulda_buffer_2_5deg'

#define some basic data
occ.asc = read.asc(occ_grid)                 # read in the vegetation grid
#rf.asc = occ.asc; rf.asc[which(is.finite(rf.asc) & rf.asc!=1)] = 0  #set all veg != 1 (rainforests) to 0

################################################################################
#get a subset of the data for occur & background
pos = as.data.frame(which(is.finite(occ.asc),arr.ind=TRUE)) #get all points that have data
pos$occ = occ.asc[cbind(pos$row,pos$col)] #append the occurrence data
pos.subset = pos[which(pos$row %% 2 == 0 & pos$col %% 2 == 0),] #get a subset of the data as a regular grid of every second cell (even row / col numbers)
#pos.subset= pos  #this line is an alternative to the previous, which selects every 2nd cell
pos.subset = rbind(pos.subset,pos[which(pos$occ==1),])
pos.subset = unique(pos.subset) #ensure all occurrence cells are included in dataset

#append the current environmental data
for (tfile in list.files(current.bioclim,pattern='\\.asc.gz',full.name=TRUE)) {
	tasc = read.asc.gz(tfile) #read in the data
	dataname = gsub(current.bioclim,'',tfile); dataname = gsub('\\.asc.gz','',dataname); dataname = gsub('/','',dataname) #define the column name
  pos.subset[dataname] = tasc[cbind(pos.subset$row,pos.subset$col)] #append the data
}
pos.subset = na.omit(pos.subset) #ensure there is no missing data
# one site (of 28) is lost that is outside the landmask but has all the env data.

#now remove the mask column (if any)
if (exists("mask_layer_name")) {
  pos.subset[mask_layer_name] <- NULL
}

#define the occurrences & background ... then write out the data
occur = data.frame(species=species_name,pos.subset[which(pos.subset$occ==1),])
occur$occ = NULL #define the occurrences
bkgd = data.frame(species='bkgd',pos.subset); bkgd$occ=NULL #define the background

write.csv(occur,'occur.csv',row.names=FALSE) #write out the occurrences
write.csv(bkgd,'bkgd.csv',row.names=FALSE) #write out the background

################################################################################
#run maxent
dir.create(output_folder_name)
output_folder_name <- paste(work.dir,output_folder_name,sep="")
threads_term = paste("threads=",threads,sep="")
sys_command = paste('java -mx2048m -jar ',maxent.jar,' -e bkgd.csv -s occur.csv -o', output_folder_name,'nothreshold nowarnings novisible -P -J -r -a ',threads_term)
system(sys_command)
