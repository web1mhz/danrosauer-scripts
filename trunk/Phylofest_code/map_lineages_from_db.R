# this script will extract records from the Moritz Lab database and generate maps and kml files
rm(list=ls())

library(RODBC)
library(maptools)
library(sp)
library(ggmap)

# parameters
genus <- 'Heteronotia'
species <- 'binoei'
buffer <- 3
aus_xmax <- 154
aus_xmin <- 112
aus_ymax <- -9
aus_ymin <- -44
output_prefix <- 'H_binoei_'
output_dir <- "C:/Users/u3579238/Dropbox/ARC Laureate/maps/lineage maps/Heteronotia"
#################

ch <- odbcConnect(dsn="Moritz_specimen_ANSI")
odbcGetInfo(ch)
specimen_columns <- sqlColumns(ch,"specimen",special=F)

sqlQuery <- paste("select * from genus where genus = '",genus,"'",sep="")

genus_row  <- sqlQuery(ch,query=sqlQuery)
genus_id <- genus_row$genus_id
rm(genus_row)

sqlQuery <- paste("select * from specimen where genus_id = ",genus_id," and species = '",species,"' and lineage_from_sequence is not NULL ",
                  "and lineage_from_sequence != '' and latitude is not NULL and longitude is not NULL",sep="")
records <- sqlQuery(ch,query=sqlQuery,stringsAsFactors=F)
odbcClose(channel=ch)

lineages <- unique(records$lineage_from_sequence)

proj <- "+proj=longlat +datum=WGS84"
kml_icon = "http://maps.gstatic.com/mapfiles/ridefinder-images/mm_20_red.png"

for (lineage in lineages) {
  
  lineage_records <- records[records$lineage_from_sequence==lineage,]
  lineage_XY <- lineage_records[,c("longitude","latitude")]
  lineage_record_count <- nrow(lineage_records)
  cat("\nLineage:",lineage,"\t\trecords:",lineage_record_count,"\n")
  
  # make a map
  #windows(10,10)
  filename <- paste(output_dir,"/",output_prefix,lineage,".png",sep="")
  png(filename,width=1024,height=768)
  xmin <- max(min(lineage_records$longitude) - buffer,aus_xmin)
  xmax <- min(max(lineage_records$longitude) + buffer,aus_xmax)
  ymin <- max(min(lineage_records$latitude) - buffer,aus_ymin)
  ymax <- min(max(lineage_records$latitude) + buffer,aus_ymax)
  
  #map_url <- get_map(source="google",location=c(xmin,ymin,xmax,ymax),maptype="terrain",urlonly=TRUE,api_key="test")
  #map_url <- paste(map_url,"key=AIzaSyBlhkwZIDowVSPZLZ2eG_T2841R0GvD8T4",sep="")
  map <- get_map(source="google",location=c(xmin,ymin,xmax,ymax),maptype="terrain")
  
  m <- ggmap(map,extent='panel',main=lineage)
  m <- m + geom_point(aes(x=longitude,y=latitude), data = lineage_XY, size = 4, alpha=0.5,color="red")
  m <- m + ggtitle(lineage)
  m <- m + theme(axis.title=element_text(face="bold", size="12"))
  m <- m + theme(plot.title=element_text(face="bold", size="16"))
  m <- m + xlab("longitude") + ylab("latitude") 

  print(m)
  
  dev.off()

  # write a kml
  filename.kml <- paste(output_dir,"/",output_prefix,lineage,".kml",sep="")
  sp_points <- SpatialPointsDataFrame(coords=lineage_XY,data=lineage_records)
  kmlPoints(sp_points,kmlfile=filename.kml,kmlname=lineage,name=lineage_records$tissue_number, icon=kml_icon)

  cat("Done",filename.kml,"\n")

}

graphics.off()