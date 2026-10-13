bond <- read.table('bond.dat')
names(bond) <- c('ingot','metal','pressure')
library(nimble)
ingotCode <- nimbleCode({
  for (i in 1:N){
    pressure[i] ~ dnorm(mu[i],sd=serr)
    mu[i] <- metal[metn[i]] + ing[ingot[i]]
  }
  for (i in 1:p){
    metal[i] ~ dnorm(70,sd=100)
  }
  for (i in 1:q){
    ing[i] ~ dnorm(0,sd=sing)
  }
  serr ~ dgamma(1.1,.1)
  sing ~ dgamma(1.1,.1)
})
metn <- as.numeric(as.factor(bond$metal))
N <- length(metn)
p <- 3
q <- 7
ingotConsts <- list(p=p,q=q,N=N,ingot=bond$ingot,
                    metn=metn)
ingotData <- list(pressure=bond$pressure)
ingot.out <- nimbleMCMC(code=ingotCode,
                         constants = ingotConsts,
                         data=ingotData,
                         nchains=5,niter=202000,
                         nburnin=2000,thin=100,
                         samplesAsCodaMCMC = TRUE,
                         summary = TRUE,WAIC = TRUE,
                         monitors=c('metal','ing','serr','sing'))
ingot.out$summary$all.chains
ingot.out$WAIC
library(coda)
gelman.diag(ingot.out$samples)
allsamps <- as.matrix(ingot.out$samples)
raftery.diag(allsamps)
effectiveSize(allsamps)
allsamps[1,]
iron <- allsamps[,9]
plot(density(iron))
