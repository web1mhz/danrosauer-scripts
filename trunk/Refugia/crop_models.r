
library(SDMTools)
library(raster)

#define directories
input.dir   = 'C:/Users/u3579238/Work/Refugia/DistModels_Reside/Reptiles/'; setwd(input.dir)
output.dir  = 'C:/Users/u3579238/Work/Refugia/DistModels_Reside/Reptiles_cropped/'; setwd(output.dir)

input_dirs= list.dirs(path = input.dir, full.names = FALSE, recursive = FALSE)
output_files= list.files(path = output.dir, pattern = "*.asc.gz", recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)
in_folder = list.files(path=input_dirs[1], pattern = "1990.asc.gz", recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

first_model.asc = read.asc.gz(paste(input_dirs[1],"//",in_folder[1],sep=''))
first_model.ras = raster.from.asc(first_model.asc)
new_extent =      extent(first_model.ras)
new_extent@xmin = 140

for (dir in input_dirs) {
  in_folder = list.files(path=dir, pattern = "1990.asc.gz", recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)  
  tfile = in_folder[1]
  filepath=paste(dir,"/",tfile,sep='')
  pathname = strsplit(dir,'//')
  dataname = unlist(pathname)[2]  #define the column name
  outname = paste(output.dir,dataname,sep="")
  
  if (! paste(dataname,".asc.gz",sep='') %in% output_files) {
    grid.asc = read.asc.gz(filepath) #read in the data
    grid.ras = raster.from.asc(grid.asc)
    grid_crop.ras = crop(grid.ras,new_extent,snap='in')
    grid_crop.asc= asc.from.raster(grid_crop.ras)
    write.asc.gz(grid_crop.asc,outname)
    cat("\nCropped asc written for",dataname)
  } else {
    cat("\nSkipped",dataname)
  }
}
