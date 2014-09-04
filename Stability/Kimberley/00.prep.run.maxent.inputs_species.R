# adapted from code by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
# this version by Dan Rosauer is designed to model using specuies locations
# queried directly from the Moritz Lab database


################################################################################
rm (list=ls())
library(SDMTools)
library(RODBC)
library(raster)

################################################################################
################################################################################

#define directories
work.dir            <- 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/'; setwd(work.dir)
current.bioclim     <- 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000'
maxent.jar          <- 'C:/Users/u3579238/Work/Refugia/Stability/maxent.jar'

genus = 'Carlia'
species = 'johnstonei'

use_only_with_lineage_ID <- FALSE

mask_layer_name     <- 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/landmask.asc'
output_folder_name  <- paste('maxent.output',genus,species,sep="_")
output_raster_name   <- paste(genus,"_",species,".asc",sep="")

maxent_threads      <- 8

ymin <- -24
resolution <- 2.5/60

####################################
# get the database location records
####################################

ch <- odbcConnect(dsn="Moritz_db_dotcloud")
odbcGetInfo(ch)
specimen_columns <- sqlColumns(ch,"specimen",special=F)

sqlQuery <- paste("select * from genus where genus = '",genus,"'",sep="")

genus_row  <- sqlQuery(ch,query=sqlQuery)
genus_id <- genus_row$genus_id

sqlQuery <- paste("select latitude, longitude from specimen where genus_id = ",genus_id," and species = '",species,
                  "' and (is_observation = 0  or is_observation is NULL) and latitude is not NULL",
                  " and longitude is not NULL ", "and (usable_spatial = 1 or usable_spatial is NULL)",sep="")
if (use_only_with_lineage_ID) {
  sqlQuery <- paste(sqlQuery, " and (lineage_from_mtDNA <> '' and lineage_from_mtDNA is not NULL)")
}

records <- sqlQuery(ch,query=sqlQuery,stringsAsFactors=F)
odbcClose(channel=ch)

coords <- unique(records[,c("longitude","latitude")])
coords <- round(coords/resolution) * resolution
mask.ras <- raster(mask_layer_name)

species.ras <- mask.ras
species.ras[species.ras[] > 0] <- 0
cells <- unique(cellFromXY(species.ras, coords))
species.ras[cells] <- 1
#writeRaster(species.ras,output_raster_name)

species.asc <- asc.from.raster(species.ras)
rm(species.ras)

cat("\n\n*** Starting maxent model with", nrow(coords),"unique locations of", genus, species,"***\n\n")

################################################################################
#get a subset of the data for occur & background
pos = as.data.frame(which(is.finite(species.asc),arr.ind=TRUE)) #get all points that have data
pos$species = species.asc[cbind(pos$row,pos$col)] #append the species data
pos.subset = pos[which(pos$row %% 3 == 0 & pos$col %% 3 == 0),] #get a subset of the data as a regular grid of every second cell (even row / col numbers)
#pos.subset= pos  #this line is an alternative to the previous, which selects every 2nd cell
pos.subset = rbind(pos.subset,pos[which(pos$species==1),]); pos.subset = unique(pos.subset) #ensure all cells in the target area are included in dataset

# apply a southern limit
 xy <- getXYcoords(species.asc)
 ycoords <- which(xy$y >= ymin)
 pos.subset <- pos.subset[pos.subset$col %in% ycoords,]

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
  #pos.subset["landmask"] <- NULL
}

#define the occurrences & background ... then write out the data
#occur = data.frame(species='species',pos.subset[which(pos.subset$species==1),])
occur           <- pos.subset[which(pos.subset$species==1),]
occur$species   <- NULL
species_col     <- rep(paste(genus,"_",species,sep=""),nrow(occur))
occur <- cbind(species=species_col, occur)
rm(species_col)

bkgd = data.frame(species='bkgd',pos.subset)
bkgd$species=NULL #define the background

write.csv(occur,'occur.csv',row.names=FALSE) #write out the occurrences
write.csv(bkgd,'bkgd.csv',row.names=FALSE) #write out the background

################################################################################
#run maxent
if (! length(dir(output_folder_name)) > 0) {  # create the maxent output folder if it doesn't already exist
  dir.create(output_folder_name)
}

sys_command = paste('java -mx2048m -jar ',maxent.jar,' -e bkgd.csv -s occur.csv -o', output_folder_name,'nothreshold nowarnings novisible -P jackknife -r -a') # variant with jacknifing for better info on variable importance

# to run in parallel
if (maxent_threads > 1) {
  sys_command <- paste(sys_command, " threads=",maxent_threads,sep="")
}

system(sys_command)
