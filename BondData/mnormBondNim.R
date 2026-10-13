library(nimble)
mnormBondCode <- nimbleCode({
  for (i in 1:nIng){
    y[i,1:3] ~ dmnorm(mu[i,1:3],cov=tau[1:3,1:3])
    for (j in 1:nTmt){
      mu[i,j] <- b0[trt[i,j]]
    }
  }
  for (i in 1:nTmt){
    b0[i] ~ dnorm(70, sd=10)
  }
  tau[1:3,1:3] ~ dinvwish(s[1:3,1:3],3)
  
})
nIng=7
nTmt=3
dat <- read.table('bond.dat')
names(dat) <- c('ingot','metal','strength') 
datord <- dat[order(dat$ingot),]  # order by ingot
# put in matrix shape
dim(datord)
#  number of ingots is 21/3
resp <- t(matrix(datord$strength,3,7)) # row is response for each ingot.
treatment <- t(matrix(datord$metal,3,7))  
s=matrix(6,3,3)
diag(s) <- 12
s <- chol(s)
mnormConsts <- list(nIng=nIng,nTmt=nTmt,
                    trt=treatment,
                    s=s)
mnormData <- list(y=resp)
mnormBond.out <- nimbleMCMC(code=mnormBondCode,
                        data=mnormData,
                        constants=mnormConsts,
                        nchains=8,niter=102000,
                        nburnin=2000,thin=100,
                        samplesAsCodaMCMC = TRUE,
                        summary=TRUE,WAIC=TRUE,
                        monitors=c('b0','tau'))
mnormBond.out$summary
mnormBond.out$WAIC
library(coda)
gelman.diag(mnormBond.out$samples[,1:3])
allsamps <- as.matrix(mnormBond.out$samples)
raftery.diag(allsamps)
plot(density(allsamps[,1]-allsamps[,2]))
plot(density(allsamps[,1]-allsamps[,3]))
plot(density(allsamps[,2]-allsamps[,3]))
mean(allsamps[,1]-allsamps[,2] <0)
mean(allsamps[,1]-allsamps[,3] <0)
mean(allsamps[,2]-allsamps[,3] <0)
