# create group_species_lineage_.csv

rm(list=ls())

base.dir  <-  'C:/Users/u3579238/Work/AMT/Models/'
work.dir  <-  paste(base.dir,"lineage_sites",sep="")
genus     <-  'Heteronotia'
#in.file.name    <-  paste(genus,"_lin_loc.csv",sep="")
in.file.name <- "Heteronotia_lin_loc_from_db_18mar14_nth_of_22S.csv"
out.file.name   <-  paste(base.dir,"Heteronotia_binoei_lineage_list.csv",sep="")

setwd(work.dir)

tab <- read.csv(in.file.name)
tab <- tab[tab$use==1,c("Genus","Species","ModelGroup","lineage")]
tab <- unique(tab)
lin_num <- seq(from=1,to=nrow(tab),1)
tab <- cbind(tab,lin_num)

write.csv(tab,out.file.name,row.names=F)

# now number lineages in lin_loc
lin_loc <- read.csv(in.file.name)
lin_loc <- lin_loc[lin_loc$use == 1,]
lin_loc_merged <- merge(lin_loc,tab[,c("lineage","lin_num")],by='lineage')
#lin_loc_merged <- lin_loc_merged[,c(2:5,1,6:11,14)]
#lin_loc_merged <- lin_loc_merged[,c(2:5,1,6:17)]
#write.csv(lin_loc_merged,paste(genus,"_lin_loc2.csv",sep=""),row.names=FALSE)
write.csv(lin_loc_merged,"Heteronotia_lin_loc_from_db_18mar14_nth_of_22S_numbered.csv",row.names=FALSE)
