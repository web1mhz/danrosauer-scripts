######################METRICS
#require(entropart)


##all metrics combine

##AED/HEAD functions without ecoPD
#from (http://rfunctions.blogspot.com/2013/09/functions-for-phylogenetic-diversity.html)
#spp<-colnames(samp)
get.nodes1<- function(tree, spp){
	edge<- which.edge(tree, spp)
	nodes<- tree$edge[edge,1] 
	root.edge<- which(tree$edge[,1]==(length(tree$tip.label)+1))
while(!(edge %in% root.edge)){
	edge<- which.edge(tree, tree$edge[edge,1])
	next.node<- tree$edge[edge,1]
	nodes<- c(nodes, next.node)
	}
	nodes
	}

#tree and comm data must match
AED<-function(abund,tree){
	
p.tree<-multi2di(tree, random = TRUE)
uninodes<-p.tree$edge[,2][-which(p.tree$edge[,2] %in% 1:length(p.tree$tip.label))]
infonode<-cbind(p.tree$edge,p.tree$edge.length)
nodevalues<-numeric(length(uninodes))
for(i in 1:length(uninodes)){
temptree<-extract.clade(p.tree,uninodes[i])
abundvalues<-abund[which(names(abund) %in% temptree$tip.label)]
nodevalues[i]<-(infonode[which(infonode[,2]==uninodes[i]),3])/(sum(abundvalues))
}

nodeabval<-cbind(uninodes,nodevalues)
tipnums<-data.frame(cbind(p.tree$tip.label,1:length(p.tree$tip.label)))
AED<-numeric(length(abund))
for(j in 1:length(abund)){
sp<-tipnums$X1[j]
nodes <- get.nodes1(p.tree,as.vector(sp))
splength<-infonode[which(infonode[,2]==as.numeric(as.vector(tipnums[which(tipnums[,1]==sp),2]))),3]
splength<-splength/(abund[which(names(abund)==sp)])
AED[j]<-(sum(nodeabval[which(nodeabval[,1] %in% nodes),2])) + splength
}
names(AED)<-tipnums[,1]
AED<-unlist(AED)
return(AED)
}
#AED(samp,tree)



#for PD functions, IAC, IC, gamma, HED, HAED, Chaos
#IC and Gamma standardized for comparison across dif size trees
pd_fcts<-function(samp, tree){
	
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute pd")
    }
    
    species <- colnames(samp)
    SR <- rowSums(ifelse(samp > 0, 1, 0))
    nlocations = dim(samp)[1]
    nspecies = dim(samp)[2]
    PDs = NULL
    PDavgs=NULL
    IACs=NULL
    HEDs=NULL
    HAEDs=NULL
    ICs=NULL
    GAMs=NULL
    Chao0s=NULL
    Chao1s=NULL
    Chao2s=NULL
    PENTs=NULL
    ScheinerDPs=NULL
    ScheinerDPabs=NULL
    
    for (i in 1:nlocations) {
        present <- species[samp[i, ] > 0]
        treeabsent <- tree$tip.label[which(!(tree$tip.label %in%present))]
        subsite<-samp[i,which(samp[i,]>0)]
        #if no species present    
        if (length(present) == 0) {
            PD <- 0
            PDavg<-0
            IAC<-NA
            HED<-NA
            HAED<-NA
            IC<-NA
            GAM<-NA
            Chao0=NA
            Chao1=NA
            Chao2=NA
            PENT=NA
            ScheinerDP=NA
            ScheinerDPab=NA

            }else{
        	
        	if (length(present) == 1) {
        	 #sum of single present branch
                PD <- node.age(tree)$ages[which(tree$edge[, 2] ==which(tree$tip.label == present))]
                PDavg<-PD/1 
                IAC<-NA  
                HED<-NA
                HAED<-NA
                IC<-NA
                GAM<-NA
                Chao0=NA
                Chao1=NA
                Chao2=NA
                PENT=NA
                ScheinerDP=NA
                ScheinerDPab=NA
              }else{

       if (length(treeabsent) == 0) {
            #all branches present
            PD <- sum(tree$edge.length)
            PDavg<-PD/SR
            IAC<-sum(sapply(subsite,function(x){abs(x-(sum(subsite)/length(subsite)))}))/tree$Nnode
            HED<--1*sum((evol.distinct(tree,"fair.proportion")$w/PD)*log(evol.distinct(tree,"fair.proportion")$w/PD))
           HAED<- -1*sum(((subsite[match(evol.distinct(tree,"fair.proportion")$Species,names(subsite))]*AED(subsite,tree))/PD)*log((subsite[match(evol.distinct(tree,"fair.proportion")$Species,names(subsite))]*AED(subsite,tree)/PD)))
            
            IC<-foo(tree,metric="colless")
            GAM<-foo(tree,metric="gamma")
            relabund<-subsite/sum(subsite)
            Chao0=ChaoPD(as.numeric(relabund),0, as.hclust.phylo(tree), Normalize=TRUE)
            Chao1=ChaoPD(as.numeric(relabund),1, as.hclust.phylo(tree), Normalize=TRUE)
            Chao2=ChaoPD(as.numeric(relabund),2, as.hclust.phylo(tree), Normalize=TRUE)
            PENT=PhyloEntropy(relabund,1, as.hclust.phylo(tree), Normalize=TRUE)$Total
            
            ScheinerDP=sum((evol.distinct(tree,"fair.proportion")$w/PD)^2)^(1/(1-2))
            ScheinerDPab=sum((relabund*evol.distinct(tree,"fair.proportion")$w/PD)^2)^(1/(1-2))

                    
              }else{
            sub.tree <- drop.tip(tree, treeabsent)
            sub.tree.depth <- max(node.age(sub.tree)$ages)
            orig.tree.depth <- max(node.age(tree)$ages[which(tree$edge[,2] %in% which(tree$tip.label %in% present))])
            
            PD <- sum(sub.tree$edge.length) + (orig.tree.depth - sub.tree.depth)
            PDavg<-PD/length(present)      
            IAC<-sum(sapply(subsite,function(x){abs(x-(mean(subsite)))}))/sub.tree$Nnode
            HED<--1*sum((evol.distinct(sub.tree,"fair.proportion")$w/PD)*log(evol.distinct(sub.tree,"fair.proportion")$w/PD))
           GAM<-foo(sub.tree,metric="gamma")
            
            relabund<-subsite/sum(subsite)
 
            ScheinerDP=sum((evol.distinct(sub.tree,"fair.proportion")$w/PD)^2)^(1/(1-2))
            ScheinerDPab=sum((relabund*evol.distinct(sub.tree,"fair.proportion")$w/PD)^2)^(1/(1-2))
			
   			if(length(subsite)==2){
   				IC<-NA
   				PENT<-NA
   				Chao0<-NA
   				Chao1<-NA
   				Chao2<-NA
   				HAED<-NA
   				
   			}else{
   			IC<-foo(sub.tree,metric="colless")
   			HAED<- -1*sum(((subsite[match(evol.distinct(sub.tree,"fair.proportion")$Species,names(subsite))]*AED(subsite,sub.tree))/PD)*log((subsite[match(evol.distinct(sub.tree,"fair.proportion")$Species,names(subsite))]*AED(subsite,sub.tree)/PD)))
   			#match community  data and tree data
   			#relabund1<-as.numeric(relabund)
   			Chao0=ChaoPD(relabund,0,as.hclust.phylo(sub.tree), Normalize=TRUE)
            		Chao1=ChaoPD(relabund,1,as.hclust.phylo(sub.tree), Normalize=TRUE)
            		Chao2=ChaoPD(relabund,2,as.hclust.phylo(sub.tree), Normalize=TRUE)
            		#names(relabund1)<-names(subsite)
			PENT=PhyloEntropy(relabund,1,as.hclust.phylo(sub.tree),Normalize=TRUE)$Total
            }
          }
	}
}
        
PDs <- c(PDs, PD)
PDavgs<-c(PDavgs,PDavg)
IACs<-c(IACs,IAC)
HEDs<-c(HEDs,HED)
HAEDs<-c(HAEDs,HAED)
ICs<-c(ICs,IC)
GAMs<-c(GAMs,GAM)
Chao0s=c(Chao0s,Chao0)
Chao1s=c(Chao1s,Chao1)
Chao2s=c(Chao2s,Chao2)
PENTs=c(PENTs,PENT)
ScheinerDPs=c(ScheinerDPs,ScheinerDP)
ScheinerDPabs=c(ScheinerDPabs,ScheinerDPab)
    }
    
    PDout <- data.frame(PD = PDs, PDavg = PDavgs, SR = SR, IAC=IACs,HED=HEDs,HAED=HAEDs,IC=ICs,GAM=GAMs,Chao0=Chao0s,Chao1=Chao1s,Chao2=Chao2s,PENT=PENTs,ScheinerDP=ScheinerDPs,ScheinerDPab=ScheinerDPabs)
    rownames(PDout) <- rownames(samp)
    return(PDout)
	}
#pd_fcts(samp,tree)

###Endemism functions
#require(phylobase)

parent.prob <- function(probabilities) {
  # probabilities is a vector of values between 0 and 1
  # add code to check values are of correct type!
  parent.prob <- 1 - prod(1-probabilities)
  return(parent.prob)
}

scale.to <- function(vec,vec.sum) {
  #mat is a vector
  #this function rescales each the vector values to sum to 'sum'
  vec.tot <- sum(vec,na.rm=TRUE)
  if (vec.tot > 0) {
    vec.out <- vec.sum*vec/vec.tot
  } else {
    vec.out <- rep(0,times=length(vec))  #should columns for tips / branches with no occurrences be removed?
  }
  return(vec.out)
}

## ---------------
##For PE and PD(abund)
endemism_fcts<-function (samp, tree) {

### step 1: trimming the tree to match the community data table thus creating a "community tree" (sensu Cam Webb).

if (length(tree$tip.label) > ncol(samp)) {
	tree <- drop.tip (tree, which(!tree$tip.label %in% colnames(samp))) }

### step 2: converting a community tree into a MRP matrix

phylomatrix <- matrix (0, length(tree$tip.label), length(tree$edge.length))
for (i in 1:length(tree$tip.label)) {
	lineage <- which (tree$edge[,2] == i)
	node <- tree$edge[lineage,1]
	while (node > length(tree$tip.label)+1) {
		branch <- which (tree$edge[,2] == node)
		lineage <- c (lineage, branch)
		node <- tree$edge[branch,1]
		}
	phylomatrix[i,lineage] = 1
	}

# this script fills a matrix of 0's with 1's corresponding to the incidence of an OTU on a particular branch.
# the code is pretty slow on large trees.

### step 3: re-ordering the OTUs to match.

phylomatrix <- phylomatrix[sort.list(tree$tip.label), ]
samp <- samp[ ,sort.list(colnames(samp))]
samp <- as.matrix(samp)

### step 4: creating a community phylogeny matrix indicating incidence of branches per site

commphylo <- samp %*% phylomatrix
commphylo <- ifelse (commphylo > 0, 1, 0)

commphyloAb <- samp %*% phylomatrix


### step 5: calculate the range (in sites) of each branch

ranges <- colSums(commphylo)
rangesAb <- colSums(commphyloAb)

### step 6: calculate weightings per branch

weightsPE <- t(commphylo) / ranges
weightsPE[which(ranges==0),]<-0
weightsPE <- ifelse(weightsPE>0,1,0)

weightsAb<-t(commphylo) / rangesAb
weightsAb[which(rangesAb==0),]<-0


### step 6: calculating PD per sample
SR <- rowSums(ifelse(samp > 0, 1, 0))
pe <- t(weightsPE) %*% tree$edge.length
pd.ab<- t(weightsAb) %*% tree$edge.length
pdab_avg<-pd.ab/SR
output<-cbind(pe,pd.ab,pdab_avg)
colnames(output)<-c("PE","PD_Ab","PD_Ab_Avg")
return(output)

}

#######
#MPD functions

weighted.var <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
na.rm)
}

mpd_fcts<-function(samp,tree){
	
    N <- dim(samp)[1]
    mpd <- numeric(N)
    mpdAb<-numeric(N)
    vpd<-numeric(N)
    vpdAb<-numeric(N)
    mntdAb<-numeric(N)
    vntdAb<-numeric(N)
    mntd<-numeric(N)
    vntd<-numeric(N)
    
    dis<-cophenetic(tree)
    
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
  			
  			#mean dist
            sample.weights <- t(as.matrix(samp[i, sppInSample, 
                  drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                  drop = FALSE])
             mpdAb[i] <- weighted.mean(sample.dis, sample.weights)
          	 vpdAb[i] <- weighted.var(sample.dis,sample.weights)
             mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
             vpd[i] <- var(sample.dis[lower.tri(sample.dis)])
             
             #min dist
             diag(sample.dis)<-NA
             mntds <- apply(sample.dis, 2, min, na.rm = TRUE)
             sample.weights2 <- samp[i, sppInSample]
             mntdAb[i] <- weighted.mean(mntds, sample.weights2)
             vntdAb[i] <- weighted.var(mntds, sample.weights2)
             mntd[i]<-mean(apply(sample.dis,2,min,na.rm=TRUE))
             vntd[i]<-var(apply(sample.dis,2,min,na.rm=TRUE))

        }
        else {
            mpd[i] <- NA
            mpdAb[i] <- NA
            vpd[i] <- NA
            vpdAb[i] <- NA
            mntdAb[i] <- NA
            vntdAb[i] <- NA
            mntd[i]<-NA
            vntd[i]<-NA

        }
    }
    return(cbind(MPD=mpd,MPD_Ab=mpdAb,VPD=vpd,VPD_Ab=vpdAb,MNTD_Ab=mntdAb,VNTD_Ab=vntdAb,MNTD=mntd,VNTD=vntd))
    
}

rao_fcts<-function(samp,tree){
	
	raoent<-raoD(samp,tree)$Dkk #equal to MPDab
	simp<-raoD(samp)$Dkk
	intraMPDab<-raoent*simp
	return(cbind(Raos=raoent,intraMPDab=intraMPDab))
	
}
