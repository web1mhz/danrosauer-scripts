
library(SDMTools)
library(raster)

#define directories
input.dir   = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models/'; setwd(input.dir)
output.dir  = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models_ext/'
template_ext ='C:/Users/u3579238/GISData/EnvironmentGrids/AusGDMGrids/ForMaxent/bio1.asc'

input_files= list.files(path = input.dir, full.names = FALSE, recursive = FALSE)
output_files= list.files(path = output.dir, pattern = "*.asc.gz", recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

template.asc = read.asc(template_ext)
template.ras = raster.from.asc(template.asc)
new_extent =      extent(template.ras)
new_extent@xmin = max(new_extent@xmin, 140.01)

for (tfile in input_files) {
  checkname=unlist(strsplit(tfile,".",fixed=T))
  if (checkname[length(checkname)]=="asc") {   # only accept filenames ending in .asc
    filepath=paste(input.dir,tfile,sep='')
    pathname = tfile
    dataname = unlist(strsplit(pathname,".asc")[1])[1]  #define the column name
    outname = paste(output.dir,dataname,sep="")
    
    if (! paste(dataname,".asc",sep='') %in% output_files) {
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
  } else {
    cat("\nSkipped",tfile)
  }
}
