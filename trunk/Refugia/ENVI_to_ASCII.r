library(raster)
library(rgdal)

setwd("C:\\Users\\u3579238\\GISData\\Mackey\\fpar_stats\\totfgreen_coefvar")
ras <- raster("tot_fgreen_coefvar_042000-042012")
plot(ras)

writeRaster(ras, "tot_fgreen_coefvar_042000-042012.asc", "ascii")

setwd("C:\\Users\\u3579238\\GISData\\Mackey\\fpar_stats\\mean_stdev_fgreen")
ras2 <- raster("tot_fgreen_mean_stdev_042000-042012")
writeRaster(ras, "tot_fgreen_mean_stdev_042000-042012.asc", "ascii")