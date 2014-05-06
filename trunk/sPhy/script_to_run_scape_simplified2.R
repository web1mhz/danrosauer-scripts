#Nov. 22, 2013. Matt Helmus, Caroline Tucker, Arne Mooers, Lanna Jin
##Script to run "SCAPE" modules: analysis of phylogenetic metrics using simulated landscapes and trees.

##The simplified version omits the perturbations and missing information modules. Will add to later versions...

##Load Libraries
require(plyr)
require(apTreeshape)
require(phytools)
require(ape)
require(geiger)
require(phangorn)
require(picante)
require(entropart)

##Necessary source files
setwd("C:/Users/u3579238/Work/Software/danrosauer-scripts/sPhy")
source("tree_generate_fct.R")
source("gradient_cart_simulation_v3.1.R")
source("landscape_manipulation_plotting.R")
source("combined_metrics1.1.R")


#functions to combine lists for results file
f<-function(x){function(i){sapply(x,`[[`,i)}}
f2<-function(x){as.data.frame(Map(f(x),names(x[[1]])))}


#############################FUNCTIONS

##MODULE 1 - generate landscapes
runscape<-function(numsp=64,scape.type=c("all",1,2,3,4,5,6,7,8),scape.size=15,treelist,maxabund=100,plotscapes=FALSE)
{

	#Scape types (second number is example # on figure, which is probably confusing): 
	
	##Environment covary w phylogeny
	#1 = #1a) ##lambda E = 1, lambda R = 0, repulsion
	#weak attraction (i.e. env filtering), strong repulsion(competition), with no signal of range size

	#2 = #1b) ##lambda E = 1, lambda R = 0, attraction
	#strong attraction (i.e. env filtering), weak repulsion (competition), with no signal of range size

	##Environment and range covary w phylogeny
	#3 = #2a)##lambda E = 1, lambda R = 1, repulsion
	
	#4 = #2b)##lambda E = 1, lambda R = 1, attraction
	
	##No environment or range covariance with phylogeny
	#5 = #3 spatial)##lambda E =0, lambda R = 0
	
	#6 = #3 aspatial)##lambda E =0, lambda R = 0
	
	##Range phylogeny covariance with no environmental
	#signal
	#7 = #4 spatial)  lambda E =0, lambda R = 1
	
	#8 = #4 aspatial)##lambda E =0, lambda R = 1,
	  
  ##############Generate objects to hold & identify data###########
  
  #number of trees input
  ntrees<-length(treelist)
  
  #types of landscapes called for, see Figure
  if(any(scape.type=="all"))
  {scape.type=c(1:8)}
  scapenames<-c("ex1a","ex1b","ex2a","ex2b","ex3spatial","ex3aspatial","ex4spatial","ex4aspatial")
  scapenames<-scapenames[scape.type]
  
  
  #Parameters for landscapes, lists
  g.center<-list(0.2,20,0.2,20,1,1,1,1)
  g.range<-list(1,1,0.2,20,1,1,1,1)
  g.repulse<-list(0.2,1,0.2,1,1,1,1)
  repulse<-list(TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE)
  signal.center<-list(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)
  signal.range<-list(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,TRUE)
  same.range<-list(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE)
  
  Spmatrixlist=list(NULL)	
  ###########Start simulating landscapes########################
  
  #create landscapes and metrics for all input trees	
  for(m in 1:ntrees){
  	
  scapelist<-lapply(scape.type,function(x){scape(treelist[[m]],scape.size=scape.size,g.center=g.center[[x]],g.range=g.range[[x]],g.repulse=g.repulse[[x]],repulse=repulse[[x]],wd.all=30,signal.center=signal.center[[x]],signal.range=signal.range[[x]],same.range=same.range[[x]],site.stoch.scale=1,sd.center=1,sd.range=1.5,th=3)})	
  names(scapelist)<-scapenames			
  		
  species<-colnames(scapelist[[1]]$Y)
  
  ################Prepare landscapes into siteXspecies matrices, calculate abundances for abundance weightings, plot 2-D landscape if desired####
  
  spmatrixlist<-rep(list(NA),length(scape.type))
  names(spmatrixlist)<-paste("mat",scapenames,sep="_")
  
  spmatrixlist<-lapply(1:length(scape.type),function(x){c(landscape_abunds(scapelist[[x]],maxabund,species,figs=plotscapes),TREE=m)})
  Spmatrixlist<-c(Spmatrixlist,spmatrixlist)
  
  }
  
  names(Spmatrixlist)<-paste("mat",scapenames,sep="_")
  
  return(Spmatrixlist[-1])
	
}

##OUTPUT is treelist and list of landscapes##
##END MODULE 1

###############Apply perturbations#######################
#None in this simplified version



############
##MODULE 2##
#####################Calculate metrics######################

scapemetrics<-function(Spmatrixlist,treelist){
	
  species<-colnames(Spmatrixlist[[1]]$abundance.matrix)	

  ##calculate metrics from combined_metrics script
  rescompleteinfo <-lapply(Spmatrixlist,FUN=function(y){
  	colnames(y$abundance.matrix)<-species
  	tr<-y$TREE
  	a<-pd_fcts(y$abundance.matrix,treelist[[tr]])
  	b<-endemism_fcts(y$abundance.matrix,treelist[[tr]])
  	c<-mpd_fcts(y$abundance.matrix,treelist[[tr]])
  	d<-rao_fcts(y$abundance.matrix,treelist[[tr]])
  	cbind(a,b,c,d,tr)
	})

  ########output metric results
  
  results<-do.call(rbind.data.frame,rescompleteinfo)
  ident<-rownames(results)
  identspl<-do.call(rbind.data.frame,(strsplit(ident,"[.]")))
  colnames(identspl)<-c("Landscape","CellID")
  results2<-cbind(identspl,results)
  return(results2)
}

###############USING SCAPE SETUP
####################Generate trees
numsp=64
numtrees=4
tree_gen<-tree_generate(numsp,numtrees)
treeinfo<-tree_gen$treeinfo
treelist<-tree_gen$treelist

###################Generate scapes
#default function values
scape.size=15 
total.cells=(scape.size+1)^2
maxabund=100
subsample.percent=75
detect.threshold=5
plotscapes=FALSE
scape.type="all"

Spmatrixlist<-runscape(numsp=64,scape.type="all",scape.size=15,treelist,maxabund=100,plotscapes=FALSE)


#################Calculate metrics

output<-scapemetrics(Spmatrixlist,treelist)
write.table(output,file="output.txt")
a<-matrix(NA,nrow=27,ncol=1)

for(i in unique(output$Landscape)){
  a<-cbind(a,cor(subset(output,output$Landscape==i)[,3:29],use="complete")[,1])
}

a<-a[,-1]
apply(a,1,"var")