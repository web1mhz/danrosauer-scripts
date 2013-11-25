# create group_species_lineage_.csv

base.dir  <-  'C:/Users/Dan/Work/AMT/Models/'
work.dir  <-  paste(base.dir,"lineage_sites",sep="")
genus     <-  'Gehyra'
in.file.name    <-  paste(genus,"_lin_loc.csv",sep="")
out.file.name   <-  paste(base.dir,"group_lineage_list.csv",sep="")

setwd(work.dir)

tab <- read.csv(in.file.name)
tab <- tab[,c("Genus","Species","ModelGroup","lineage")]
tab <- unique(tab)

write.csv(tab,out.file.name,row.names=F)

