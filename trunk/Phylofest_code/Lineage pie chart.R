
library(raster)
library(sp)
library(stringr)

#define directories
input.dir       <- 'C:/Users/u3579238/Work/AMT/Models/lineage_models/asc_clipped'
output.dir      <- 'C:/Users/u3579238/Work/AMT/Maps/lineage_models/Heteronotia/'
file_pattern    <- '[:alnum:]*Heteronotia_b.*asc$'
group_lin_file  <- 'C:/Users/u3579238/Work/AMT/Models/group_lineage_list.csv'
preface         <- 'lin_model_'
suffix          <- ''

lineage_site_file   <- 'C:/Users/u3579238/Work/AMT/Models/lineage_sites/Heteronotia_lin_loc_from_db_18mar14_nth_of_22S.csv'

too.low = 0.05

####  end of parameters

files = list.files(path = input.dir, pattern = file_pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

lin_sites <- read.csv(lineage_site_file)
lin_sites <- subset(lin_sites,long > 0)

setwd(input.dir)

lin.stack <- stack(files)

i <- 1
lineage <- vector(mode="character")

splitfiles <- strsplit(files,"_")
for (file in splitfiles) {
  this.lineage <- file[5]
  lineage[i] <- strsplit(this.lineage,"[.]")[[1]][1]
  i <- i+1
}
rm(i)

names(lin.stack) <- lineage

x <- 134
y <- -16.4
xy <- SpatialPoints(data.frame(x,y))

lin_spot_values <- extract(lin.stack,xy)

lin_spot_values <- lin_spot_values / sum(lin_spot_values) # make them sum to 1
names(lin_spot_values) <- lineage

# collapse the small ones to 'other'
values <- as.numeric(lin_spot_values)
minors <- which(values < too.low)
main.values <- lin_spot_values[,-minors]
n <- length(main.values)
main.values[n+1] <- sum(lin_spot_values[,minors])
names(main.values)[n+1] <- paste("Others <",(round(too.low,3) * 100),"%",sep="")
header <- paste("lineages at ",x,y)

windows(12,6)
par(mfrow=c(1,2))
pie(main.values, radius=0.9, col=rainbow(n+1),labels=names(main.values),main=header)

header <- paste("lineages at ",x,y,"without <",too.low)
pie(main.values[1:n], radius=0.9, col=rainbow(n+1),labels=names(main.values[1:n]),main=header)
