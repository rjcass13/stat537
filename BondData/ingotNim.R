bond <- read.table('bond.dat')
names(bond) <- c('ingot','metal','pressure')
library(nimble)
ingotCode <- nimbleCode({
  for (i in 1:N){
    pressure[i] ~ dnorm(mu[i],sd=serr)
    mu[i] <- metal[metn[i]] + ing[ingot[i]]
  }
  # Random effects have priors
  for (i in 1:p){
    metal[i] ~ dnorm(70,sd=100)
  }
  # Random effects have hyperpriors
  # Most important thing for random effect is to get the right variance
  # We don't want that variance to be fixed, allow it to explore it
  for (i in 1:q){
    ing[i] ~ dnorm(0,sd=sing)
  }
  serr ~ dgamma(1.1,.1)
  sing ~ dgamma(1.1,.1)
})

curve(dgamma(x, 1, 1), from = 0, to = 50)

# Define Parameters
metn <- as.numeric(as.factor(bond$metal))
N <- length(metn)
p <- 3
q <- 7
ingotConsts <- list(p=p,q=q,N=N,ingot=bond$ingot,
                    metn=metn)
ingotData <- list(pressure=bond$pressure)

# Run the Nimble code
ingot.out <- nimbleMCMC(code=ingotCode,
                         constants = ingotConsts,
                         data=ingotData,
                         nchains=5,niter=202000,
                         nburnin=2000,thin=100,
                         samplesAsCodaMCMC = TRUE,
                         summary = TRUE,WAIC = TRUE,
                         monitors=c('metal','ing','serr','sing'))

# Examine the output
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
