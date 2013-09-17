rm(list=ls())
library(raster)

#define directories
input.dir   = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models/'; setwd(input.dir)
output.dir  = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models_aligned/'
template_ext ='C:/Users/u3579238/GISData/EnvironmentGrids/AusGDMGrids/ForMaxent/bio1.asc'

file.pattern       <- '.+asc$'  #regex
input_files= list.files(path = input.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE,
                        ignore.case = TRUE, include.dirs = FALSE)
output_files= list.files(path = output.dir, pattern = file.pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

template.ras = raster(template_ext)
new_extent =      extent(template.ras)
new_extent@xmin = max(new_extent@xmin, 140.01)

new_only <- FALSE  # if true, skip grids which are already in the output directory

for (tfile in input_files) {
#for (tfile in input_files[85:length(input_files)]) {
  filepath=paste(input.dir,tfile,sep='')
  outname = paste(output.dir,tfile,sep="")

  if ((!tfile %in% output_files) | (! new_only)) {
    grid.ras = raster(tfile)
    grid_ext.ras = extend(grid.ras,new_extent,value=0) # extend to the union of current grid and new extent
    grid_ext.ras = crop(grid_ext.ras,new_extent) # crop back to new extent
    writeRaster(grid_ext.ras,outname,overwrite=TRUE, NAflag=-9999)
    cat("\nExtended asc written for",tfile)
  } else {
    cat("\nSkipped",tfile)
  }
}
