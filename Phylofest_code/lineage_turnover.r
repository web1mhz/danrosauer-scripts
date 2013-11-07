rm(list=ls())

library(ecodist)

models <- read.csv("C:/Users/u3579238/Work/Refugia/Results/outputs_25Sep/rf_lin_models.csv")

#strip out empty cells
models <- models[which(models$total > 0),]

#group cells
group_by <- 0.1
models$lat_group <- as.integer(models$lat/group_by) * group_by
models$long_group <- as.integer(models$long/group_by) * group_by
models_only <- cbind(models$lat_group,models$long_group,paste(models$lat_group,"_",models$long_group,sep=""),models[,12:112])
names(models_only)[1:3] <- c("lat_group","long_group","cell_name")

models_grouped <- aggregate(models_only[,4:ncol(models_only)],by=list(models_only$lat_group,models_only$long_group,models_only$cell_name),FUN=sum)
names(models_grouped)[1:3] <- c("lat_group","long_group","cell_name")

# set NA values to 0
for (col in 1:ncol(models_grouped)) {
  na_count <- length(which(is.na(models_grouped[,col])))
  cat(names(models_grouped)[col], "\t",na_count,"\n")
  if (na_count>0) {
    models_grouped[which(is.na(models_grouped[,col])),col] <- 0
  }
}

row.names(models_grouped) <- models_grouped$cell_name

dist_matrix <- distance(as.matrix(models_grouped[,4:ncol(models_grouped)]), method="bray-curtis",spweight="absence")
dist_matrix <- as.matrix(dist_matrix)
gc()

dist_matrix[which(is.na(dist_matrix))] <- 1

my.nmds <- nmds(as.dist(dist_matrix),maxdim=3)
my.axes <- as.matrix(my.nmds[[1]])

#rescale to (roughly) 0 to 1
my.axes_r <- my.axes
my.axes_r[which(is.na(my.axes_r))] <- -1

rescale <- function(in_val, old_min, old_max, new_min = 0, new_max = 1) {
  out_val <- inv_val + old_min + new_min
  out_val <- (out_val / (old_max - old_min)) * new_max
  return(out_val)
}

vec.r <- unlist(my.axes[,1][1])
vec_range <- range(vec.r)
vec.r2 <- sapply(X=vec.r, FUN=rescale, old_min = vec_range[1],old_max = vec_range[2])

rowcount <- length(vec.r)
vec.r <- vec.r + rep(1,rowcount)
vec.r <- vec.r / rep(2,rowcount)
vec.g <- unlist(my.axes[,1][2])
vec.g <- vec.g + rep(1,rowcount)
vec.g <- vec.g / rep(2,rowcount)
vec.b <- unlist(my.axes[,1][3])
vec.b <- vec.b + rep(1,rowcount)
vec.b <- vec.b / rep(2,rowcount)

gc()
cat(range(vec.r),range(vec.g),range(vec.b))
my.cols = rgb(red=vec.r,vec.g,vec.b)

windows(10,10)
plot(models_grouped$long_group,models_grouped$lat_group,col=my.cols)