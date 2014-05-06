#Caroline Tucker
#Determine abundances, plot landscape

#original landscape = output from scape function

landscape_abunds<-function(original_landscape,abundmax,species,figs=TRUE)

{	
	require(vegan)
	
#calculate abundances
abundmax<-abundmax

#has one more row than the original scape demand, fix?
PA_mat<-as.matrix(original_landscape$Y)

site.size=nrow(PA_mat)

abund_mat<-original_landscape$Y*original_landscape$X.joint
abund_mat2<-as.matrix(abund_mat*abundmax)


#make latlong distances
latlong<-original_landscape$index
#cbind(rep(1:sqrt(site.size),each=sqrt(site.size)),seq(1:sqrt(site.size)))
colnames(latlong)<-c("XDim","YDim")

##Descriptive statistics && Plotting

#print mantel test results
#PA
#if(site.size<1000){
#spatial<-mantel(dist(latlong),vegdist(PA_mat,"jaccard"))
#Abundance
#mantel(dist(latlong),vegdist(PA_mat,"bray"))
#}
#else{spatial=NULL}

if(figs){
	par(ask=TRUE)
#B&W
#bw_richnessmap<-image(matrix(rowSums(PA_mat),sqrt(site.size),sqrt(site.size),byrow=TRUE),col=gray((5:0)/8),xaxt="n",yaxt="n",main="Species Richness")

col_richnessmap<-heatcol<-(colorRampPalette(c("yellow","red")))
image(matrix(rowSums(PA_mat),sqrt(site.size),sqrt(site.size),byrow=TRUE),col=heatcol(15),xaxt="n",yaxt="n",main="Species Richness")

#species x area
par(mfrow=c(1,2))
hist(colSums(PA_mat),ylab="Number of species",xlab="Number of sites",main="Species Area Relationship",col="lightgrey")
hist(colSums(abund_mat),ylab="Number of species",xlab="Number of individuals",main="Species Abundance Relationship",col="lightgrey")

#Beta diversity for random cell across range
#if(site.size<1000){
#betadist<-as.matrix(vegdist(PA_mat,"jaccard"))
#comp_beta<-as.matrix(betadist[,1])

#par(mfrow=c(1,1))
#beta_map<-image(matrix(comp_beta,sqrt(site.size),sqrt(site.size),byrow=TRUE),col=heatcol(25),xaxt="n",yaxt="n",main="Jaccard disimmilarity between cell 1 (1,1) and all other cells")#possibly not interesting
#}	


	}

colnames(abund_mat2)<-species
colnames(PA_mat)<-species
return(list(abundance.matrix=abund_mat2,PA_mat=PA_mat))

}




###Apply perturbations

#landscape.abund.matrix - abundance.matrix output from landscape_abunds

landscape_perturb<-function(landscape.abund.matrix,subsample.percent,detect.threshold){

original<-landscape.abund.matrix
nsites<-nrow(original)


if(subsample.percent<100){

original[sample(nsites,nsites*(subsample.percent)/100),]<-NA	

}

perturbed<-original

if(detect.threshold!="all"){
	perturbed[perturbed<detect.threshold]<-0
}


return(list(original_abund=original,perturbed_abund=perturbed))	
	
}

##convert abundance.matrix as necessary to PA for calculations
siteXsp_prep <- function(siteXsp, plot_names_in_col1=TRUE){

	#Move plot names in column 1 (if specified as being present) into the 
	#row names of the matrix and drop the column of names
	if(plot_names_in_col1){
		row.names(siteXsp) <- siteXsp[,1]
		siteXsp <- siteXsp[,-1]
	}
	
	#Remove any redundant species columns (i.e. species that were never observed).
	colsums <- colSums(siteXsp)
	siteXsp <- siteXsp[,colsums > 0]

	#Make the siteXsp matrix into a pres/abs. (overwrites intitial siteXsp 
	#matrix):
	siteXsp <- ifelse(siteXsp > 0, 1, 0)

	return(siteXsp)
}


#gens<-sapply(1:8,FUN=function(x){Ancestors(tree,x,"parent")})
#gens<-cbind(tree$tip.label,as.numeric(gens))
#trial<-rbind(test1a$abundance.matrix,as.numeric(gens[,2]))

#Make species level branches of phylogeny have no information = reduced info for sampling

treetogenus<-function(tree){
	tree_lowinfo<-tree
	nspecies<-length(tree$tip)
	gens<-sapply(1:nspecies,FUN=function(x){which(tree$edge[,2]==x)})
	tree_lowinfo$edge.length[gens]<-0.001
	return(tree_lowinfo)
}



#perturb trees

##record outputs

#tree stats, landscape stats, metric values, columns for each perturbation