library(ape)
library(plotrix)

#PARAMETERS
#tree = phylo object
#scape.size = edge dimension of square landscape
#g.center = phylogenetic signal in species optimal values  phylogenetic attraction (range centers)   See corBlomberg in ape, 1=brownian,<1=rates of evol accelerate, >1=rates decelerate.
#g.range = phylogenetic signal in species niche widths  (range sizes). 1=brownian,<1=rates of evol accelerate, >1=rates decelerate.
#g.repulse = include phylogenetic repulsion
#wd.all = niche width denominator, larger values make larger on average range sizes
#signal.center = T/F simulate with phylosignal in range centers
#signal.range = T/F simulate with phylosignal in range sizes
#same.range = T/F make all range sizes equal
#repulse = T/F include phylogenetic repulsion
#center.scale = adjust the strength of phylogenetic attraction independent of signal (generally leave = 1 this is more for error checking)
#range.scale = adjust the strength of phylogenetic signal in niche width (range) (generally leave = 1 this is more for error checking)
#repulse.scale = adjust the strength of phylogenetic repulsion (generally leave = 1 this is more for error checking)
#site.stoch.scale = adjust the strength of random variation in "carrying capacity" across sites
#sd.center = sd in rnorm() for the range centers, increase to get more variation in center values across species
#sd.range =  sd in rnorm() for the range sizes, increase to get more variation in range sizes across gradients
#rho = vary the overall strength of signal in the phylogeny (Grafen branch adjustment) (probably more for error checking, so leave as NULL)
#th = probability threshold 10^-th above which species are considered present at a site (increase this value to obtain more species at sites in Y (p/a matrix) and thus increase average site SR) 

scape<-function(tree, scape.size=10, g.center=1, g.range=1, g.repulse=1, wd.all=150, signal.center=TRUE, signal.range=TRUE, same.range=TRUE,repulse=TRUE,center.scale = 1, range.scale = 1, repulse.scale = 1, site.stoch.scale = .5, sd.center=1, sd.range=1,rho=NULL, th=8)
{
  #deal with the tree
    if (is(tree)[1] == "phylo")
    {
        if (is.null(tree$edge.length))
        {
            tree <- compute.brlen(tree, 1)    #Note this assigns arbitrary branch-lengths
        }
        V <- vcv.phylo(tree, corr = TRUE)     #Note this will mess with trees that are not ultrametric, such as those with argitrary branch-lengths
    } else {
        V <- tree
    }
    Vinit<-V       
  #initialize
    nspp <- dim(V)[1]
    bspp2 <- NULL
    Vcomp <- NULL

        Xscale <- 1          #scale the strength of the probability matrix X
        Mscale <- site.stoch.scale        #scale stochasticity in niche distributions 
        Vscale1 <- center.scale         #scale the strength of the optimal values on axis one
        Vscale2 <- center.scale         #scale the strength of the optimal values on axis two
    
    # Grafen's rho adjust strength of phylogenetic signal overall.
    if(!is.null(rho)){
     V <- 1-(1-Vinit)^rho
     V <- V/max(V)
    }
 
    #########################################################################################################################
    #SIMULATION
    nsites<-scape.size                                                #number of sites for the square landscape
    mx <- t(as.matrix((-(nsites)/2):(nsites/2)))                      #gradient
    #mx<-t(as.vector(scale(1:nsites)))
    m <- length(mx)                                                   #new number of sites
      
    ############
    #ATTRACTION (RANGE CENTERS,NICHE OPTIMA) AND RANGE WIDTHS
    if(signal.center){
      g<-abs(g.center)
      V.a<-vcv(corBlomberg(g, tree),corr=T)           #adjust phylogenetic signal
      iD <- t(chol(V.a))
    } else {
      V.a<-V
      iD <- diag(nspp)
    }
    if(signal.range){
      g<-abs(g.range)
      V.w<-vcv(corBlomberg(g, tree),corr=T)           #adjust phylogenetic signal
      iD.w <- t(chol(V.w))
    } else {
      V.w<-V
      iD.w <- diag(nspp)
    }
  ##environmental/geographical gradient 1   
  #if (envirogradflag1 == 1)
  #{
    e <- iD %*% rnorm(nspp,sd=sd.center)                                                               #assign optimal values 
    #e <- iD %*% runif(nspp)                                                               #assign optimal values 
    #e <- iD %*% rep(1,nspp)                                                               #assign optimal values 
    e <- Vscale1 * (e - mean(e))/apply(e, 2, sd)                                           #z-scores and optional scaling of the optimal values
    bspp1 <- e
    if(!same.range){
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          #wd <- range.scale*iD.w %*% rnorm(nspp,sd=sd.range)
          wd <- range.scale*iD.w %*% rnorm(nspp,sd=sd.range)
          #wd <- range.scale*iD.w %*% runif(nspp)
          #wd<-wd.all*wd
          wd<-wd+(abs(min(wd)))
          wd<-wd/max(wd)
          #wd[wd==0]<-sort(wd)[2]-mean(sort(wd)[-1]-sort(wd)[-length(wd)])             #this might cause errors, removed the zero in wd and replace it with nonzero minimum minus the average gap between the sorted wd values.
          #wd[wd==0]<-sort(wd)[2]                                                       #Assign the zero to the nonzero minimum
          dif<-sort(wd)[-1]-sort(wd)[-length(wd)]
          rat<-mean(dif/sort(wd)[-1])
          wd[wd==0]<-sort(wd)[2]-sort(wd)[2]*rat                                #Assign the zero with the mean rato of nearist neighbor distances over the larger item
          wd<-wd.all*wd
          X <- exp(-((spmx - mxsp)^2)/t(matrix(wd,nspp,m))) #Niche distributions
        } else {
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
        #X <- 1/(1 + exp(-(b0scale * array(1, c(m, 1)) %*% rnorm(nspp) + mx %*% t(e))))        #alternative calculation
          X <- exp(-((spmx - mxsp)^2)/wd.all) #Niche distributions
    }       
    X <- Xscale * X                                                                       #Scales this initial species x site probability matrix 
    #Xsmooth <- X                                                                         #Distributions without random variation
    X1 <- diag(1 - Mscale * runif(m)) %*% X                                               #Scale and include random variation into the niche distributions
    #}
    
    ##environmental/geographical gradient 2
    #if (envirogradflag2 == 1) {
        e <- iD %*% rnorm(nspp,sd=sd.center)
        e <- Vscale2 * (e - mean(e))/apply(e, 2, sd)
        bspp2 <- e
        #mx2 <- as.matrix(mx[sample(m)])                                   #make second gradient random
        #mx2 <- t(as.matrix((-(nsites)/2):(nsites/2)))                      #gradient construct
        #mx2<-t(as.matrix(mx)) 
     if(!same.range){
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          wd <- range.scale*iD.w %*% rnorm(nspp,sd=sd.range)
          wd<-wd+(abs(min(wd)))
          wd<-wd/max(wd)
          #wd[wd==0]<-sort(wd)[2]-mean(sort(wd)[-1]-sort(wd)[-length(wd)])             #this might cause errors, removed the zero in wd and replace it with nonzero minimum minus the average gap between the sorted wd values.
          #wd[wd==0]<-sort(wd)[2]                                                       #Assign the zero to the nonzero minimum
          dif<-sort(wd)[-1]-sort(wd)[-length(wd)]
          rat<-mean(dif/sort(wd)[-1])
          wd[wd==0]<-sort(wd)[2]-sort(wd)[2]*rat                                #Assign the zero with the mean rato of nearist neighbor distances over the larger item
          wd<-wd.all*wd
          #wd <- iD.w %*% rnorm(nspp,sd=sd.range)
          #wd <- range.scale*(wd-mean(wd))/sd(as.vector(wd))
          X <- exp(-((spmx - mxsp)^2)/t(matrix(wd,nspp,m))) #Niche distributions     
     } else {
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
         #X <- 1/(1 + exp(-(b0scale * array(1, c(m, 1)) %*% rnorm(nspp) + mx %*% t(e))))        #alternative calculation
          X <- exp(-((spmx - mxsp)^2)/wd.all) #Niche distributions
     }
     X <- Xscale * X
     #Xsmooth <- Xsmooth * X
     X2 <- diag(1 - Mscale * runif(m)) %*% X
     #}

     ##################################        
     #REPULSION
     X.repulse <- NULL
      if (repulse) {
        compscale <- repulse.scale
        b0scale <- 0
        g<-abs(g.repulse)
        V.r<-vcv(corBlomberg(g, tree),corr=T)          #adjust phylogenetic signal 
      #calculate the repulsion matrix
        Vcomp <- solve(V.r, diag(nspp))
        Vcomp <- Vcomp/max(Vcomp)
        Vcomp <- compscale * Vcomp
        iDcomp <- t(chol(Vcomp))
        colnames(Vcomp) <- rownames(Vcomp)
        bcomp <- NULL
        for (i in 1:m) {
          bcomp <- cbind(bcomp, iDcomp %*% rnorm(nspp))
        }
        bcomp0 <- 0
        Xcomp <- exp(bcomp0 + bcomp)/(1 + exp(bcomp0 + bcomp))
        #X <- X * t(Xcomp)
        X1<-X1 * t(Xcomp)
        X2<-X2 * t(Xcomp)
        X.repulse<-t(Xcomp)
      }
      
##################################
#JOINT PROBABILITY MATRIX
  X.<-NULL
  spp.Xs<-array(NA,dim=c(m,m,nspp))
  for(i in 1:nspp){
    sppX<-matrix((X1[,i]) %*% t(X2[,i]))
    spp.Xs[,,i]<-sppX
    X.<-cbind(X.,matrix(sppX))
   }
   colnames(X.)<-colnames(X2)
   #X.attract <- X
  ######################
  #PA matrix
  m.<- dim(X.)[1]
  Y <- matrix(0, ncol = nspp, nrow = m.)
  #Y[matrix(runif(nspp * m.), ncol = nspp) < X.] <- 1
  Y[10^-th < X.] <- 1                                        #could also use a hard threshold
  colnames(Y) <- colnames(X.)
  index<-NULL
  #m<-length(mx)
  index<-cbind(matrix(sapply(1:m,rep,times=m)),matrix(rep(1:m,times=m)))
  colnames(index)<-c("X1","X2")
  
  ########### OUTPUT  
    return(list(Y = Y, index=index, 
                X.joint=X., X1=X1, X2=X2, sppXs=spp.Xs, 
                V.phylo=Vinit, V.phylo.rho = V, V.center = V.a, V.range = V.w, V.repulse = Vcomp, 
                bspp1 = bspp1, bspp2 = bspp2, u = mx, wd=wd.all))
  
}  #function end

#VALUES
#Y = presence/absence matrix
#index = spatial coordinates for X and Y (stacked columns)
#X.joint = full probabilities of species at sites, used to construct Y 
#X1 = probabilities of species along gradient 1
#X2 = probabilities of species along gradient 2
#sppXs = full probabilities of each species as an array arranged in a scape.size X scape.size matrix  
#V.phylo = initial phylogenetic covariance matrix from tree, output of vcv.phylo(tree, corr=T)
#V.phylo.rho = phylogenetic covariance matrix from tree scaled by Grafen if rho is provided, otherwise just an output of vcv.phylo(tree, corr=T)
#V.center = scaled (by g.center) phylo covariance matrix used in the simulations
#V.range =  scaled (by g.range) phylo covariance matrix used in the simulations
#V.repulse = scaled (by g.repulse) phylo covariance matrix used in the simulations
#bspp1 = species optima for gradient 1
#bspp2 = species optima for gradient 2                               
#u = the env gradients values for the two gradients
#wd = the denominator for species ranges 

############################
#EXAMPLE AND SOME PLOTS
#tree<-stree(8,type="balanced")       #signal in centers
#kk<-scape(tree, scape.size=100, g.center=100, g.range=1, g.repulse=1, wd.all=150, signal.center=TRUE, signal.range=FALSE, same.range=FALSE, repulse=FALSE,center.scale = 1, range.scale = 1, repulse.scale = 1, site.stoch.scale = 0, sd.center=3, sd.range=1,rho=NULL, th=20)

#graphics.off()
#par(mfrow=c(1,Ntip(tree)),mar=c(.1,.1,.1,.1))
#for(j in 1:Ntip(tree)){
#color2D.matplot(1 - kk$sppXs[,,j]/max(kk$sppXs[,,j]), xlab = "", ylab = "",main = "",border=NA,do.hex=FALSE,axes=FALSE)
#}

#windows()
#par(mfrow=c(2,1))
#matplot((kk$X1), type = "l", xlab="gradient",ylab = "probability", main = "Gradient 1",col=rainbow(dim(kk$X1)[2]),lty=1)
#matplot((kk$X2), type = "l", xlab="gradient",ylab = "probability", main = "Gradient 2",col=rainbow(dim(kk$X2)[2]),lty=1)

#windows()
#plot(x=1:dim(kk$Y)[1],y = rowSums(kk$Y), main = "SR",type = "l")
#cor(kk$X1)
