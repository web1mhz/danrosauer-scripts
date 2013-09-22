
rm(list=ls())

# parameters #############################################################
input.dir   = 'C:/Users/u3579238/Work/Refugia/Results/outputs_11Jun'; setwd(input.dir)
output.dir  = 'C:/Users/u3579238/Work/Refugia/Results/outputs_11Jun'

lineage_ranges <- read.csv("lineage_ranges_0.9.csv")
region_areas   <- read.csv("region_areas.csv")
lineage_ranges_end <- lineage_ranges[lineage_ranges$endemic==TRUE,]
region_names = c("CYP","AWT","MEQ","CEC","SEA")

windows(8,8)
plot(region_areas$X,log(region_areas$region_area),col="red",ylim=c(5.5,11))
points(lineage_ranges_end$end_region,log(lineage_ranges_end$range),xlab = "Region", ylab = "Log lineage range")

windows(8,8)
boxplot(log(range)~end_region,data=lineage_ranges_end)
plot(region_areas$X,log(region_areas$region_area),col="red",ylim=c(5.5,11))
points(lineage_ranges_end$end_region,log(lineage_ranges_end$range),xlab = "Region", ylab = "Log lineage range")

lin_models_occ <- read.csv("rf_lin_models_occ_0.9.csv")
lineages <- names(lin_models_occ)[12:ncol(lin_models_occ)]

# calculate lineage stats
col_count <- ncol(lineage_ranges)
lin_count <- nrow(lineage_ranges)
col_of_zeroes <- rep(0,lin_count)
lineage_stats <- cbind(lineage_ranges,col_of_zeroes,col_of_zeroes,col_of_zeroes)
names(lineage_stats)[(col_count+1):(col_count+3)] <- c("median_stability","percentile_10_stability","conf_interval")
for (j in 1:length(lineages)) {
  lineage <- as.character(lineage_stats$lineage[j])
  cells <- which(lin_models_occ[,lineage]==1)
  lineage_stats$median_stability[j] <- median(lin_models_occ[cells,"stabil_static"],na.rm=T)
  lineage_stats$percentile_10_stability[j] <- quantile(lin_models_occ[cells,"stabil_static"],0.1,na.rm=T)
  lineage_stats$conf_interval[j] <- quantile(lin_models_occ[cells,"stabil_static"],0.975,na.rm=T) - quantile(lin_models_occ[cells,"stabil_static"],0.025,na.rm=T)
}

median_order <- order(lineage_stats$median_stability,decreasing=TRUE)

# plot in order of decreasing median stability
windows(8,8)
plot(0,0,xlim=c(1,lin_count),ylim=c(0,1),col=NULL,xlab="Lineage", ylab="Stability", main="Ordered by median stability") # set up an empty plot
i=0
for (j in median_order) {
  i <- i+1
  lineage = as.character(lineage_stats$lineage[j])
  median = lineage_stats$median_stability[j]
  perc_10 =lineage_stats$percentile_10_stability[j]
  cells <- which(lin_models_occ[,lineage]==1)
  region <- lin_models_occ[cells,"region"]
  values <- lin_models_occ[cells,"stabil_static"]
  colours <- region
  colours[which(values< perc_10)] <- grey(0.7)
  points(rep(i,length(cells)),values,pch=20,col=colours)
  points(i,median,col="purple", pch=15,cex=0.8)
  points(i,perc_10,col="yellow", pch=15,cex=0.8)  
}
legend(0.05,1.03,legend=c(region_names,"below 10th percentile"),fill=c(1:5,grey(0.7)),cex=0.8,bty='n')
points(c(13,13),c(0.97,0.92),col=c("purple","yellow"),pch=15)
text(c(14,14),c(0.97,0.92),pos=4,labels=c("Median","10th percentile"),cex=0.9)

# order by 10th percentile of stability (ie only found in stable v less picky)
perc_10_order <- order(lineage_stats$percentile_10_stability,decreasing=TRUE)

windows(8,8)
plot(0,0,xlim=c(1,lin_count),ylim=c(0,1),col=NULL,xlab="Lineage", ylab="Stability", main="Ordered by 10th percentile of stability") # set up an empty plot
i=0
for (j in perc_10_order) {
  i <- i+1
  lineage = as.character(lineage_stats$lineage[j])
  median = lineage_stats$median_stability[j]
  perc_10 =lineage_stats$percentile_10_stability[j]
  cells <- which(lin_models_occ[,lineage]==1)
  region <- lin_models_occ[cells,"region"]
  values <- lin_models_occ[cells,"stabil_static"]
  colours <- region
  colours[which(values< perc_10)] <- grey(0.7)
  points(rep(i,length(cells)),values,pch=20,col=colours)
  points(i,median,col="purple", pch=15,cex=0.8)
  points(i,perc_10,col="yellow", pch=15,cex=0.8)  
}

legend(0.05,1.03,legend=c(region_names,"below 10th percentile"),fill=c(1:5,grey(0.7)),cex=0.8,bty='n')
points(c(13,13),c(0.97,0.92),col=c("purple","yellow"),pch=15)
text(c(14,14),c(0.97,0.92),pos=4,labels=c("Median","10th percentile"),cex=0.9)

# order by 95% confidence interval of stability (ie found in a wide range of stability v not)
conf_interval_order <- order(lineage_stats$conf_interval,decreasing=TRUE)

windows(8,8)
plot(0,0,xlim=c(1,lin_count),ylim=c(0,1),col=NULL,xlab="Lineage", ylab="Stability", main="Ordered by size of 95% confidence interval") # set up an empty plot
i=0
for (j in conf_interval_order) {
  i <- i+1
  lineage = as.character(lineage_stats$lineage[j])
  median = lineage_stats$median_stability[j]
    cells <- which(lin_models_occ[,lineage]==1)
  region <- lin_models_occ[cells,"region"]
  values <- lin_models_occ[cells,"stabil_static"]
  colours <- region
  colours[which(values > quantile(values,0.975))] <- grey(0.7)
  colours[which(values < quantile(values,0.025))] <- grey(0.7)
  points(rep(i,length(cells)),values,pch=20,col=colours)
  points(i,median,col="purple", pch=15,cex=0.8)
}

legend(0.05,1.03,legend=c(region_names,"outside 95% confidence interval"),fill=c(1:5,grey(0.7)),cex=0.8,bty='n')
points(13,0.97,col="purple",pch=15)
text(14,0.97,pos=4,labels="Median",cex=0.9)

lineages <- names(lin_models_occ)[12:ncol(lin_models_occ)]
windows(8,8)
plot(0,0,xlim=c(3.2,10.1),ylim=c(0,1),col=NULL,xlab="log Range Size", ylab = "Median Stability") # set up an empty plot
for (lineage in lineages) {
  cells <- which(lin_models_occ[,lineage]==1)
  range=length(cells)
  points(rep(log(range),range),lin_models_occ[cells,"stabil_static"],pch=20)
  points(log(range),median(lin_models_occ[cells,"stabil_static"]),col="red",pch=20,cex=1.5)  
}
