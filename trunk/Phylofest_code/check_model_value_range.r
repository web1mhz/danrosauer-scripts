

library(SDMTools)

#define directories
base.dir   = 'C:/Users/u3579238/Work/Refugia/DistModels_Reside/Reptiles_cropped/'; setwd(base.dir)
#base.dir   = 'C:/Users/u3579238/Work/Phylofest/Models/combined/lineage_models'; setwd(base.dir)
output.dir = 'C:/Users/u3579238/Work/Refugia/Results/'
file_suffix = ".asc"  # frog models end in _1999.asc.gz, reptile models don't

richness_output = "rept_rich_lin" #"frog_rich_sp_thr0.5"
endemism_output = "rept_end_lin"  #"frog_end_sp_thr0.5"

threshold = 0.0  # this is not a species level threshold, but a generic threshold across all the models

files = list.files(path = base.dir, pattern = file_suffix, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

# trying a list of coords based on the grid size
#pos = data.frame(generate_coords(model_rows,model_cols))

for (tfile in files) {
  if (length(strsplit(tfile,".asc"))==1) {                    #only use files whose name ends in .asc
    tasc = read.asc.gz(tfile)                                 #read in the data
    dataname=gsub(file_suffix,'',tfile)
    cat("\n",dataname,min(tasc,na.rm=T),max(tasc,na.rm=T))
  }
}
