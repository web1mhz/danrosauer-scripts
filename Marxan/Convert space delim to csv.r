base_dir <- "C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/marxan_batch_21"

setwd(base_dir)

spec <- read.fwf("spec.dat",c(5,10,9,21))
names(spec) <- c("id","target","spf","name")
spec <- spec[-1,]
write.csv(spec,"spec.csv")

pu <- read.fwf("pu.dat",5)
pu <- as.data.frame(pu)
pu <- as.data.frame(as.integer(pu[-1,]))
names(pu) <- "pu"
write.csv(pu,"pu.csv",quote=FALSE)

puvspr2 <- read.fwf("puvspr2.dat", c(5,5,8), stringsAsFactors=FALSE)
puvspr2 <- puvspr2[-1,]
names(puvspr2) <- c("species","pu","amount")
puvspr2$species <- as.integer(puvspr2$species)
puvspr2$pu <- as.integer(puvspr2$pu)
puvspr2$amount <- as.numeric(puvspr2$amount)
names(puvspr2) <- c("species","pu","amount")

write.csv(puvspr2,"puvspr2.csv", row.names=F)
