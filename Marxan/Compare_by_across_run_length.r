# plot the test of Marax iteration numbers for global mammals - Nov 2013

runs   <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_jan_2014/test_run_length_spf.csv")
runs_3 <- runs[which(runs$SPF==3),]
runs_2 <- runs[which(runs$SPF==2),]

windows(10,10)
par(mfrow=c(2,3))

boxplot(formula=Score~Iterations,data=runs_3,xlab="iterations",ylab="objective function", cex.axis=0.7, boxwex=0.8)
boxplot(formula=Penalty~Iterations,data=runs_3,xlab="iterations",ylab="penalty", cex.axis=0.7)
boxplot(formula=Cost~Iterations,data=runs_3,xlab="iterations",ylab="cost")

boxplot(formula=Score~Iterations,data=runs_2,xlab="iterations",ylab="objective function", cex.axis=0.7, boxwex=0.8)
boxplot(formula=Penalty~Iterations,data=runs_2,xlab="iterations",ylab="penalty", cex.axis=0.7)
boxplot(formula=Cost~Iterations,data=runs_2,xlab="iterations",ylab="cost")


windows(10,10)
par(mfrow=c(3,4))

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_500k.csv")
n <- 50
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="500k iterations",xlab="proportion",freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_1m.csv")
n <- 50
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="1m iterations",xlab="proportion", freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_2.5m.csv")
n <- 50
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="2.5m iterations",xlab="proportion", freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_5m.csv")
n <- 50
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="5m iterations",xlab="proportion", freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_10m.csv")
n <- 50
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="10m iterations",xlab="proportion", freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_20m.csv")
n <- 25
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="20m iterations",xlab="proportion", freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_40m.csv")
n <- 20
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="40m iterations",xlab="proportion",xlim=c(0,1), freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_80m.csv")
n <- 10
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="80m iterations",xlab="proportion", freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_160m.csv")
n <- 10
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="160m iterations",xlab="proportion", freq=FALSE)

ssoln <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/batch_test/test_run_length/test_run_length_low_threshold/output_ssoln_320m.csv")
n <- 5
#ssoln <- ssoln[ssoln$number>0,]
ones <- 100 * length(which(ssoln$number==n)) / nrow(ssoln)
zeroes <- 100 * length(which(ssoln$number==0)) / nrow(ssoln)
cat("zeroes",zeroes,"\tin between",100-(zeroes+ones),"\tones",ones,"\n")
hist((ssoln$number/n),n+1,main="320m iterations",xlab="proportion", freq=FALSE)
