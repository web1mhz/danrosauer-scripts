# this script will extract records from the Moritz Lab database and generate maps and kml files
rm(list=ls())

library(maptools)
library(sp)
library(raster)

# parameters

plot_species_points <- T
genus <- "Carlia"
species <- "amax"
use_only_with_lineage_ID <- F

if (plot_species_points) {
  output_prefix <- paste(genus,species,sep="_")
  base_dir  <- "C:/Users/u3579238/Work/Refugia/Stability/Kimberley/"
  models_dir <- paste(base_dir,"maxent.output_",output_prefix,"/",sep="")
  output_dir <- paste(base_dir,"maps_",output_prefix,"/",sep="")

} else {

  output_prefix <- 'Carlia amax'
  models_dir <- "C:/Users/u3579238/Work/Refugia/Stability/Kimberley/maxent.output_Carlia_amax/"
  output_dir <- "C:/Users/u3579238/Work/Refugia/Stability/Kimberley/maps_Carlia_amax"
}

per_pic   <- 1

subset <- c(1, 4, 7, 13, 19, 23, 27, 32, 42, 52, 57, 62)
do_subset=T

plot_region_boundary <- F

#################

setwd(models_dir)

models <- dir(,path=models_dir,pattern = ".asc",include.dirs = T,full.names = T)
if (do_subset) {models <- models[subset]}

coast.shp <- shapefile("C:/Users/u3579238/GISData/aus_1m.shp")
coast.shp <- coast.shp[coast.shp$NAM != "PAPUA NEW GUINEA" & coast.shp$NAM != "INDONESIA",]

if (plot_region_boundary) {
  IBRA.shp <- shapefile("C:/Users/u3579238/GISData/IBRA/IBRA7_regions.shp")
  #region_list <- "Northern Kimberley"
  #region_list <- c("Northern Kimberley", "Central Kimberley", "Victoria Bonaparte")
  #region_list <- c("Arnhem Coast","Arnhem Plateau","Central Arnhem","Daly Basin","Darwin Coastal","Gulf Fall and Uplands",
  #                  "Gulf Coastal","Mount Isa Inlier","Ord Victoria Plain","Pine Creek","Sturt Plateau","Tiwi Cobourg",
  #                  "Northern Kimberley", "Central Kimberley", "Victoria Bonaparte")
  region_list <- c("Arnhem Coast","Arnhem Plateau","Central Arnhem","Daly Basin","Darwin Coastal","Gulf Fall and Uplands",
                   "Gulf Coastal","Mount Isa Inlier","Pine Creek","Tiwi Cobourg","Northern Kimberley", "Central Kimberley",
                   "Victoria Bonaparte")
  IBRA.shp    <- IBRA.shp[IBRA.shp$REG_NAME_7 %in% region_list,]
}

if (plot_species_points) {

  ####################################
  # get the database location records
  ####################################
  library(RODBC)
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

}

pic_count <- 0

for (model in models) {

  cat("mapping ",model,"\n")
  model.ras  <- raster(model)
  model_date <- substr(model,nchar(model)-6,nchar(model)-4)
  if (do_subset) {
    filename <- paste(output_dir,"/",output_prefix," selected_",model_date,"_x",per_pic,".png",sep="")
  } else {
    filename <- paste(output_dir,"/",output_prefix,"_",model_date,"_x",per_pic,".png",sep="")
  }
  # make a map
  #windows(19.2,12)

  if (pic_count %% per_pic ==0) {
    png(filename,width=1920,height=1200)
    par(mfrow=c(2,2))
  }

  main_header <- paste(output_prefix ,model_date, sep=" ")

  plot(model.ras, main = main_header, xlab="longitude", ylab="latitude", ylim=c(-24,-10), xlim = c(118,142), cex.main=1.5)
  plot(coast.shp,lwd=0.6, add=T)
  if (plot_region_boundary) {
    plot(IBRA.shp,lwd=1.5, add=T)
  }

  if(plot_species_points) {
    points(coords,pch=20,col="blue")
  }

  if (pic_count %% per_pic == per_pic - 1) {
    dev.off()
  }
  pic_count <- pic_count + 1
}

graphics.off()
