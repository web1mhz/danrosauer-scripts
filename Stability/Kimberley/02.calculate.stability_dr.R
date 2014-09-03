#drafted by Jeremy VanDerWal ( jjvanderwal@gmail.com ... www.jjvanderwal.com )
#GNU General Public License .. feel free to use / distribute ... no warranties

rm(list=ls())

################################################################################
library(SDMTools)
library(maptools)

################################################################################
#Define some directories
work.dir = 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/maxent.output_narrow_AMT2/'; setwd(work.dir)
out.dir = 'C:/Users/u3579238/Work/Refugia/Stability/Kimberley/stability.narrow_AMT2/'; dir.create(out.dir)
Oz.shape_path = 'C:/Users/u3579238/GISData/aus5m_coast.shp'

IBRA.shp = readShapePoly("C:/Users/u3579238/GISData/IBRA/IBRA7_regions.shp")
#region_list <- "Northern Kimberley"
#region_list <- c("Northern Kimberley", "Central Kimberley", "Victoria Bonaparte")
#region_list <- c("Arnhem Coast","Arnhem Plateau","Central Arnhem","Daly Basin","Darwin Coastal","Gulf Fall and Uplands",
#                  "Gulf Coastal","Mount Isa Inlier","Ord Victoria Plain","Pine Creek","Sturt Plateau","Tiwi Cobourg",
#                  "Northern Kimberley", "Central Kimberley", "Victoria Bonaparte")
region_list <- c("Arnhem Coast","Arnhem Plateau","Central Arnhem","Daly Basin","Darwin Coastal","Gulf Fall and Uplands",
                 "Gulf Coastal","Mount Isa Inlier","Pine Creek","Tiwi Cobourg","Northern Kimberley", "Central Kimberley",
                 "Victoria Bonaparte")
IBRA.shape    <- IBRA.shp[IBRA.shp$REG_NAME_7 %in% region_list,]


#list the projections
sims = list.files(,pattern='\\.asc'); sims = gsub('\\.asc','',sims)

#define the cell size
cell.size = 0.04166667

#zero offset as cost of 1/value can produce Inf
zero.offset = 1e-7

#define the dispersal cost type "no.cost", "linear", "quadratic"
disp.type = "linear"

#define some plot info
tcols = c(colorRampPalette(c('green','brown'))(14),colorRampPalette(c('brown','yellow'))(70),colorRampPalette(c('yellow','orange'))(119),colorRampPalette(c('orange','red'))(197))
legend.pnts = cbind(c(113,114.5,114.5,113),c(-44,-44,-38,-38))
Oz.shape = readShapePoly(Oz.shape_path)

# display a progressive screen plot as it runs?
plot.live = FALSE  # this is informative, but SLOWER

################################################################################
################################################################################
# define some functions
#create images
image.grid = function(tasc,tfile,zlim=NULL,cols=NULL) {
	if (is.null(cols)) { cols = colorRampPalette(c('green','brown','yellow','orange','red'))(101) }
	#plot the image
	png(tfile,height=847/100,width=1080/100,units='cm',res=300,pointsize=5)
		layout(matrix(c(rep(rep(1,8),6),1,2,2,2,2,1,1,1),nr=7,nc=8,byrow=TRUE))
		par(mar=c(0,0,0,0),cex=1)
		if (is.null(zlim)) { zlim = range(as.vector(tasc),na.rm=T) }
		image(tasc,ann=FALSE,axes=FALSE,zlim=zlim,col=cols)
		legend.gradient(cbind(c(113,114.5,114.5,113),c(-44,-44,-38,-38)),limits=round(zlim,2),title='stability',cols=cols,cex=1.2)
		plot(Oz.shape,add=T,border="black",pbg="transparent",lwd=0.6) #add the subregions
	  plot(IBRA.shape,add=T,border="red",pbg="transparent",lwd=0.6) #add the subregions
	dev.off()
}

#calculate the moving window cost
mw.cost = function(n) {
	#define the size of the moving window
	max.dist.moved = mdd * n
	max.num.cells = round(max.dist.moved / cell.size); cat('max radius in cells is',max.num.cells,'\n')
	#create the moving window
	tt = matrix(NA,nrow=max.num.cells*2+1,ncol=max.num.cells*2+1)
	#populate the distances
	for (y in 1:(max.num.cells*2+1)){
		for (x in 1:(max.num.cells*2+1)){
			tt[y,x] = sqrt((max.num.cells+1-y)^2 + (max.num.cells+1-x)^2)
		}
	}
	#remove distance values > max.num.cells
	tt[which(tt>max.num.cells)]=NA
	#define the dispersal costs
	if (disp.type=="no.cost") {
		tt[which(is.finite(tt))]=1; tt = -log(tt)
	} else if (disp.type=="linear") {
		tt = 1-tt/max.num.cells; tt = -log(tt); tt[which(is.infinite(tt))] = -log(zero.offset)
	} else if (disp.type=="quadratic") {
		tt = (1-tt/max.num.cells)^2; tt = -log(tt); tt[which(is.infinite(tt))] = -log(zero.offset)
	}	else {exit} #Error... need to define dispersal type
	return(tt)
}
################################################################################
################################################################################
#read in and store all the data
indata = NULL
for (ii in 1:length(sims)) { sim = sims[ii]; cat(sim,'\n')
	tasc = read.asc(paste(sim,'.asc',sep='')) #read in the data
	#create an array to store the data and setup the base.asc object
	if (is.null(indata)) { indata = array(data=NA,dim=c(dim(tasc),length(sims))); base.asc = tasc; base.asc = base.asc*0 } #base.asc is for current
	indata[,,ii] = tasc #store the data
}

###create the data associated with the cost of the predicted environemtnal suitability
cost.suitability = -log(indata)
#set any Infinite values to our maximim y
cost.suitability[which(is.infinite(cost.suitability))] = -log(zero.offset)
###for calculating suitability / cost... we need to account for 'offshore' information
#to deal with this, we need to grab the maximum extent of the data and set any NA within that to 0 suitability
#that way, it will have HUGE cost and thus
pos = NULL #setup object to track largest set of values
for (ii in 1:length(sims)) { cat(sims[ii],'\n'); # extract the maximum set of points
	if (is.null(pos)) {
		pos = which(is.finite(cost.suitability[,,ii]),arr.ind=TRUE)
	} else {
		if (nrow(pos) < length(which(is.finite(cost.suitability[,,ii])))) {
      pos = which(is.finite(cost.suitability[,,ii]),arr.ind=TRUE)
		}
	}
}
pos = as.data.frame(pos) #convert to a dataframe
for (ii in 1:length(sims)) { cat(sims[ii],'\n') #cycle through, change NA's associated with pos to value of -log(zero.offset) cost
	tt = cost.suitability[,,ii][cbind(pos$row,pos$col)]
  tt = which(is.na(tt)) #define the NA positions
	cost.suitability[,,ii][cbind(pos$row[tt],pos$col[tt])] = -log(zero.offset) #set those NA's to -log(zero.offset)
}
#define the output data to store outputs
outdata = cost.suitability

#set up to do a live plot
if (plot.live) {windows(10,10)}

################################################################################
################################################################################
# calculate the OLD static predictions
min.asc = mean.asc = base.asc
#sum the predictions to get the average
for (ii in 1:length(sims)){
  cat(sims[ii],'\n')
  min.asc = pmax(min.asc,cost.suitability[,,ii],na.rm=T)
  mean.asc = mean.asc + cost.suitability[,,ii]

  # a progressive plot of the mean suitability so far, if requested
  if(plot.live) {
    plot(raster.from.asc(exp(-(mean.asc/ii))),main=paste("Static stability calculated so far, from years 000 to",sims[ii]))
  }
}

mean.asc = mean.asc / length(sims) #calculate the mean value
mean.asc = exp(-mean.asc); min.asc = exp(-min.asc) #convert back to maxent values of 0,1
write.asc(min.asc,paste(out.dir,'static.min.asc',sep=''))#write out the data
write.asc(mean.asc,paste(out.dir,'static.mean.asc',sep=''))#write out the data
bins = seq(0,1,length=101); bins = cut(0.0242,bins,labels=FALSE) # get the threshold bin for cols
cols = c(rep('gray',bins),colorRampPalette(c('brown','yellow','forestgreen'))(100)[bins:100])
image.grid(min.asc,paste(out.dir,'static.min.png',sep=''),zlim=c(0,1),cols=cols) #plot the image
image.grid(mean.asc,paste(out.dir,'static.mean.png',sep=''),zlim=c(0,1),cols=cols) #plot the image

################################################################################
################################################################################
# calculate the shifting refugia

#define cols
tcol = c(colorRampPalette(c('gray30','grey80'))(51),colorRampPalette(c('grey80','yellow','red','green','blue','saddlebrown'))(50))

# calculate the static prediction
static.asc = base.asc; static.asc[,] = cost.suitability[,,1] #set the static.asc to the first cost surface
#sum the predictions to get the average
for (ii in 2:length(sims)){
  cat(sims[ii],'\n')
  static.asc = static.asc + cost.suitability[,,ii]
}
static.asc = static.asc / (length(sims)*2-1); static.asc = exp(-static.asc)
write.asc(static.asc,paste(out.dir,'static.sum.cost.asc',sep=''))#write out the data
image.grid(static.asc,paste(out.dir,'static.sum.cost.png',sep=''),zlim=c(0,1),cols=tcol) #plot the image after logging the actual data


### 10 m/ year
#define the max dispersal distance in units per year (defined by resolution of your inputs)
#e.g., if cell size is 0.002998 decimal degrees (~250m resolution) and you want 10 m dispersal distance per year
#set mdd = 10 * 0.002998 / 250
mdd.id = 10 #m/year
mdd = mdd.id * cell.size / 4000
outdata = cost.suitability #reset the output data

#set up to do a live plot
if (plot.live) {windows(10,10)}

#calculate the stability surfaces
for (ii in (length(sims)-1):1){  #cycle through each of the layers starting with the last
  cat(sims[ii],'...')
  #define the size of the moving window
	num.years = (as.numeric(sims[ii+1]) - as.numeric(sims[ii])) *1000
	mw = mw.cost(num.years)
	#workout and store in outdata
	outdata[,,ii] = cost.suitability[,,ii] + lcmw(outdata[,,ii+1],mw,round((mdd * num.years) / cell.size))

  # a progressive plot of the mean suitability so far, if requested
  if(plot.live) {
    plot(raster(exp(-(outdata[,,ii]/(ii*2-1)))),main=paste("Dymanic stability (",mdd.id,"/yr calculated so far, from years 000 to",sims[ii]))
  }
}

#create the plot & ascii output...
tasc = base.asc; tasc[,] = outdata[,,1] #get the data
tasc[which(tasc+1>1e20)] = NA #remove the Bullshite data
tasc = tasc / (length(sims)*2-1); tasc = exp(-tasc)
write.asc(tasc,paste(out.dir,'shift.',mdd.id,'.asc',sep='')) #write out the ascii grid file
image.grid(tasc,paste(out.dir,'shift.',mdd.id,'.png',sep=''),zlim=c(0,1),cols=tcol) #plot the image after logging the actual data

### 20 m/ year
#define the max dispersal distance in units per year (defined by resolution of your inputs)
#e.g., if cell size is 0.002998 decimal degrees (~250m resolution) and you want 10 m dispersal distance per year
#set mdd = 10 * 0.002998 / 250
mdd.id = 20 #m/year
mdd = mdd.id * cell.size / 4000
outdata = cost.suitability #reset the output data

#set up to do a live plot
if (plot.live) {windows(10,10)}

#calculate the stability surfaces
for (ii in (length(sims)-1):1){ cat(sims[ii],'...') #cycle through each of the layers starting with the last
	#define the size of the moving window
	num.years = (as.numeric(sims[ii+1]) - as.numeric(sims[ii])) *1000
	mw = mw.cost(num.years)
	#workout and store in outdata
	outdata[,,ii] = cost.suitability[,,ii] + lcmw(outdata[,,ii+1],mw,round((mdd * num.years) / cell.size))
  # a progressive plot of the mean suitability so far, if requested
  if(plot.live) {
    plot(raster(exp(-(outdata[,,ii]/(ii*2-1)))),main=paste("Dynamic stability (",mdd.id,"m/yr calculated so far, from years 000 to",sims[ii]))
  }
}
#create the plot & ascii output...
tasc = base.asc; tasc[,] = outdata[,,1] #get the data
tasc[which(tasc+1>1e20)] = NA #remove the Bullshite data
tasc = tasc / (length(sims)*2-1); tasc = exp(-tasc)
write.asc(tasc,paste(out.dir,'shift.',mdd.id,'.asc',sep='')) #write out the ascii grid file
image.grid(tasc,paste(out.dir,'shift.',mdd.id,'.png',sep=''),zlim=c(0,1),cols=tcol) #plot the image after logging the actual data
