#library(ape)
#library(phytools)
#library(plyr)
#library(apTreeshape)
#library(geiger)

foo <- function(x, metric = "colless") {
    if (metric == "colless") {
        xx <- as.treeshape(x)  # convert to apTreeshape format
        colless(xx, norm="yule")  # calculate colless balance
    
    } else if (metric == "gamma") {
        gammaStat(x)
    } else stop("metric should be colless or gamma")
}

tree_generate<-function(numsp,ntrees){
	
	#generate trees
	N<-ntrees*10
	trees=rmtree(N,numsp) #simulating trees
	colless<- ldply(trees, foo, metric = "colless")  # calculate metric for each tree
	gamma<-ldply(trees, foo, metric = "gamma")
	
	#ID extreme trees/calculate quantiles
	balanced.10.percent=as.vector(quantile(colless[,1], probs=seq(0,1,0.1))[2])
	imbalanced.10.percent=as.vector(quantile(colless[,1], probs=seq(0,1,0.1))[10])
	
	#Drop the trees in the middle, generate new tree list
	trees2=list()
	j=1
	for (i in 1:N){
   	if (colless[i,1] < balanced.10.percent | colless[i,1] > imbalanced.10.percent){
   trees2[[j]] <- trees[[i]]
   j=j+1} 
	}
	
	# calculate Ic metric for each tree
	class(trees2)="multiphylo"
	Ic_trees2<- ldply(trees2, foo, metric = "colless")  

	#ultrametricize, generate chronotrees2 list - ultrametric trees with only extreme colless values
	chronotrees2=list()
	for (i in 1:length(trees2)){
	tree=chronopl(trees2[[i]], lambda=0,age.min=1,iter.max=1000,eval.max=1000)
	attr(tree,c("ploglik"))<-NULL
	attr(tree,c("rates"))<-NULL
	attr(tree,c("message"))<-NULL
	chronotrees2[[i]] <- tree
	}
	class(chronotrees2)="multiphylo"

	# calculate gamma metric for each tree
	gamma_chronotrees2<- ldply(chronotrees2, foo, metric = "gamma")  

	#Increase gamma on every tree using delta<1
	lowgammatrees <- chronotrees2
	higammatrees=list()
	for (i in 1:length(lowgammatrees)){
	treetemp <- chronotrees2[[i]]
	higammatrees [[i]]<- suppressWarnings(deltaTree(treetemp, 0.1))
		}
	class(higammatrees)="multiphylo"

	# calculate gamma metric for each extreme tree
	gamma_hi<- ldply(higammatrees, foo, metric = "gamma")  
	
	#ultrametricize average trees
	chronotrees=list()
	for (i in 1:length(trees)){
	tree=chronopl(trees[[i]], lambda=0,age.min=1,iter.max=1000,eval.max=1000)
	attr(tree,c("ploglik"))<-NULL
	attr(tree,c("rates"))<-NULL
	attr(tree,c("message"))<-NULL
	chronotrees[[i]] <- tree
	}
	class(chronotrees)="multiphylo"
	
	#final tree list - subset to the desired output
	treelist<-c(lowgammatrees,higammatrees,chronotrees)
	treeinfo<-cbind(rbind(Ic_trees2,Ic_trees2,colless),rbind(gamma_chronotrees2,gamma_hi,gamma))
	subtrees<-c(sort(sample(1:length(treelist),ntrees)))
	treelist<-treelist[subtrees]
	treeinfo<-treeinfo[subtrees,]
	colnames(treeinfo)<-c("Colless","Gamma")
	output<-list(treelist,treeinfo)
	names(output)<-c("treelist","treeinfo")
	return(output)
}

# a<-tree_generate(64,10)
# par(mfrow=c(3,3))
# lapply(a$treelist[-1],"plot")