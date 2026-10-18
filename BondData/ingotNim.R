bond <- read.table('bond.dat')
names(bond) <- c('ingot','metal','pressure')
library(nimble)
ingotCode <- nimbleCode({
  for (i in 1:N){
    pressure[i] ~ dnorm(mu[i],sd=serr)
    mu[i] <- metal[metn[i]] + ing[ingot[i]]
  }
  # Fixed effects have priors
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
# He picked the 1.1 because the shape is concave, not starting high
curve(dgamma(x, 1.1, 1), from = 0, to = 50)

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
# ICC: subject variance/ (subject + error variance)? Something like that
library(coda)
gelman.diag(ingot.out$samples)
allsamps <- as.matrix(ingot.out$samples)

# Like the Dependence factors to be close to one
raftery.diag(allsamps)

# Like ESS at least around 5k
effectiveSize(allsamps)

allsamps[1,]

# Do some comparison plots
iron <- allsamps[,9]
plot(density(iron))
lines(density(allsamps[,8]), col = "blue")

# Trace plot
plot(allsamps[, 12], type = 'l')

# Can see the long tail
plot(density(allsamps[, 12]))

# Check how often one metail is great than another
# Close to .5 implies little difference, 0 or 1 implies likelya difference
mean(allsamps[,9] - allsamps[,10]>0)
mean(allsamps[,10] - allsamps[,8]>0)


