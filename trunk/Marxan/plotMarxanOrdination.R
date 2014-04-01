rm(list=ls())

library(vegan)
library(MASS)
library(rgl)

work_dir <- 'C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/ph_step_batch/runs_combined/'
setwd(work_dir)
marxan_output <- read.csv('C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/ph_step_batch/runs_combined/all_solutions_alt_format.csv')
rownames <- paste(marxan_output$tree,marxan_output$run,sep=".")
row.names(marxan_output) <- rownames

#####################################

# less trees, for testing
marxan_output <- marxan_output[1:200,]

# exclude weird tree 10
marxan_output <- marxan_output[-(11:20),]
####################################

tree_num <- marxan_output$tree
tree_count <- length(unique(tree_num))
tree_consecutive <- 1:tree_count
tree_actual <- unique(tree_num)
tree_num_consec <- tree_num

i=0
for (num in tree_actual) {
  i=i+1
  a <- which(tree_num==num)
  tree_num_consec[a] <- tree_consecutive[i]
}

cols <- rainbow(tree_count)

numcols <- ncol(marxan_output)
marxan_output <- marxan_output[,3:numcols]

soldist<-vegdist(marxan_output,distance="bray")
sol.mds<-isoMDS(soldist)

sol.ano <- anosim(soldist, tree_num_consec)
sol.ano
plot(sol.ano)

windows(10,10)

#points
plot(sol.mds$points, pch=20, xlab='', ylab='', main='NMDS of solutions',col=cols[tree_num_consec])

# #or text
# plot(sol.mds$points, type='n', xlab='', ylab='', main='NMDS of solutions')
# text(sol.mds$points, labels=row.names(marxan_output),col=cols[tree_num])


# 3d plot
sol3d.mds<-isoMDS(soldist,k=3)

# # text
# plot3d(sol3d.mds$points, xlab = 'x', ylab = 'y', zlab = 'z', type='n', theta=40, phi=30, ticktype='detailed', main='NMDS of solutions')
# text3d(sol3d.mds$points,texts=row.names(marxan_output),pretty='TRUE',col=cols[tree_num])

#points
plot3d(sol3d.mds$points, xlab = 'x', ylab = 'y', zlab = 'z',size=1,type='s', theta=40, phi=30, ticktype='detailed', main='NMDS of solutions',col=cols[tree_num_consec])

play3d(spin3d(axis=c(1,0,0), rpm=3), duration=10)
play3d(spin3d(axis=c(0,1,0), rpm=3), duration=10)
play3d(spin3d(axis=c(0,0,1), rpm=3), duration=10)

h<-hclust(soldist, method='complete')

#bmp(file='C:\\Work\\Worldwide Marxan\\Marxan\\output\\output_dendogram.bmp',width=1900,height=1060,pointsize=10)

plot(h, xlab='Solutions', ylab='Disimilarity', main='Bray-Curtis dissimilarity of solutions',col=cols[tree_num_consec])

#dev.off()

usercut<-cutree(h,k=20)

