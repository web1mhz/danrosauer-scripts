#drafted by Dan Rosauer, using elements from a script by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )

rm(list=ls())
library(raster)

############# First a function #################################################
shp2raster <- function(shp, mask.raster, label="", value, background_value = 0,transform = FALSE, proj.from = NA,
                       proj.to = NA, map = TRUE, write_to_file = FALSE) {
# This shp2raster function adapted from: http://amywhiteheadresearch.wordpress.com/2014/05/01/shp2raster/

  require(raster, rgdal)

  # use transform==TRUE if the polygon is not in the same coordinate system as
  # the output raster, setting proj.from & proj.to to the appropriate
  # projections
  if (transform == TRUE) {
    proj4string(shp) <- proj.from
    shp <- spTransform(shp, proj.to)
  }

  # convert the shapefile to a raster based on a standardised background
  # raster
  r <- rasterize(shp, mask.raster)
  # set the cells associated with the shapfile to the specified value
  r[!is.na(r)] <- value
  # merge the new raster with the mask raster and export to the working
  # directory as a tif file
  mask.raster[!is.na(mask.raster)] <- background_value
  if (write_to_file) {
    r <- mask(merge(r, mask.raster), mask.raster, filename = label, format = "GTiff",
              overwrite = T)
  } else {
    r <- mask(merge(r, mask.raster), mask.raster)
  }

  # plot map of new raster
  if (map == TRUE) {
    plot(r, main = label, axes = F, box = F)
  }

  names(r) <- label
  return(r)
}

################################################################################

#define directories
template = 'C:/Users/u3579238/Work/Refugia/Stability/OZ.climates/bioclim/000/bioclim_01.asc'
template.ras=raster(template)

IBRA_shapefile = 'C:/Users/u3579238/GISData/IBRA/IBRA7_regions.shp'
IBRA.shp       <- shapefile(IBRA_shapefile)

# North Kimberley
region_list <- "Northern Kimberley"
IBRA_Kimb.shp  <- IBRA.shp[IBRA.shp$REG_NAME_7 %in% region_list,]
IBRA_Kimb.ras <- shp2raster(shp = IBRA_Kimb.shp, mask.raster = template.ras, value = 1, label="IBRA_kimb")

output.extent <- extent(template.ras)

IBRA_Kimb.ras <- crop(IBRA_Kimb.ras,output.extent)

output_path     <- 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/'
output_filename <- "IBRA_NorthKimberley.asc"
writeRaster(IBRA_Kimb.ras, paste(output_path,output_filename,sep=""), overwrite=TRUE)

# Broad Kimberley
region_list <- c("Northern Kimberley", "Central Kimberley", "Victoria Bonaparte")
IBRA_Kimb.shp  <- IBRA.shp[IBRA.shp$REG_NAME_7 %in% region_list,]
IBRA_Kimb.ras <- shp2raster(shp = IBRA_Kimb.shp, mask.raster = template.ras, value = 1, label="IBRA_kimb")

output.extent <- extent(template.ras)

IBRA_Kimb.ras <- crop(IBRA_Kimb.ras,output.extent)
plot(IBRA_Kimb.ras)

output_path     <- 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/'
output_filename <- "IBRA_BroadKimberley.asc"
writeRaster(IBRA_Kimb.ras, paste(output_path,output_filename,sep=""), overwrite=TRUE)

# Narrow AMT
region_list <- c("Arnhem Coast","Arnhem Plateau","Central Arnhem","Daly Basin","Darwin Coastal","Gulf Fall and Uplands",
                 "Gulf Coastal","Mount Isa Inlier","Ord Victoria Plain","Pine Creek","Sturt Plateau","Tiwi Cobourg",
                 "Northern Kimberley", "Central Kimberley", "Victoria Bonaparte")
IBRA_Kimb.shp  <- IBRA.shp[IBRA.shp$REG_NAME_7 %in% region_list,]
IBRA_Kimb.ras <- shp2raster(shp = IBRA_Kimb.shp, mask.raster = template.ras, value = 1, label="IBRA_kimb")

output.extent <- extent(template.ras)

IBRA_Kimb.ras <- crop(IBRA_Kimb.ras,output.extent)
plot(IBRA_Kimb.ras)

output_path     <- 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/'
output_filename <- "IBRA_NarrowAMT.asc"
writeRaster(IBRA_Kimb.ras, paste(output_path,output_filename,sep=""), overwrite=TRUE)
