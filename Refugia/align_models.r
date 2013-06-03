
library(SDMTools)
library(raster)

#define directories
input.dir   = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models/'; setwd(input.dir)
output.dir  = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models_aligned/'
template_ext ='C:/Users/u3579238/GISData/EnvironmentGrids/AusGDMGrids/ForMaxent/bio1.asc'

file.pattern       <- '.+asc$'  #regex
input_files= list.files(path = input.dir, pattern = file.pattern, full.names = FALSE, recursive = FALSE,
                        ignore.case = TRUE, include.dirs = FALSE)
output_files= list.files(path = output.dir, pattern = file.pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

template.asc = read.asc(template_ext)
template.ras = raster.from.asc(template.asc)
new_extent =      extent(template.ras)
new_extent@xmin = max(new_extent@xmin, 140.01)

new_only <- TRUE  # if true, skip grids which are already in the output directory

for (tfile in input_files) {
  filepath=paste(input.dir,tfile,sep='')
  pathname = tfile
  dataname = unlist(strsplit(pathname,".asc")[1])[1]  #define the column name
  outname = paste(output.dir,dataname,sep="")

  if ((!paste(dataname,".asc",sep='') %in% output_files) | (! new_only)) {
    grid.asc = read.asc(filepath) #read in the data
    grid.ras = raster.from.asc(grid.asc)
    grid_ext.ras = extend(grid.ras,new_extent,value=0) # extend to the union of current grid and new extent
    grid_ext.ras = crop(grid_ext.ras,new_extent) # crop back to new extent
    grid_ext.asc= asc.from.raster(grid_ext.ras)
    write.asc(grid_ext.asc,outname)
    cat("\nExtended asc written for",dataname)
  } else {
    cat("\nSkipped",dataname)
  }
}
