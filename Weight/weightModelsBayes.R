weight <- read.table('weight.dat',header=TRUE)
sub <- paste(substr(weight$tmt,1,1),sprintf("%02d",weight$subject),sep='')
weight$sub <- sub
ztime <- weight$time - 1
weight$ztime <- ztime
weight$tmtnum <- as.numeric(as.factor(weight$tmt))
weight$subnum <- as.numeric(as.factor(weight$sub))
library(nimble)

#########################
### Random Intercepts ###
#########################

ranintCode <- nimbleCode({
  for (i in 1:N){
    response[i] ~ dnorm(mu[i],sd=serr)
    mu[i] <- beta0[tmtnum[i]] + beta1[tmtnum[i]]*ztime[i] + u[subnum[i]] 
  }
  for (i in 1:3){
    beta0[i] ~ dnorm(80,sd=5)
    beta1[i] ~ dnorm(0,sd=1)
  }
  for (i in 1:57){
    u[i] ~ dnorm(0,sd=sdsub)
  }
  serr ~ dgamma(1.1,1)
  sdsub ~ dgamma(1.1,.5)
})
N <- 399
tmtnum <- weight$tmtnum
subnum <- weight$subnum
response <- weight$response
ztime <- weight$ztime
ranintConsts <- list(N=N,tmtnum=tmtnum,subnum=subnum,ztime=ztime)
ranintData <- list(response=response)
ranint.out <- nimbleMCMC(
  code=ranintCode,
  constants = ranintConsts,
  data=ranintData,
  nchains=5,niter=201000,
  nburnin=1000,thin=50,
  samplesAsCodaMCMC = TRUE,
  summary=TRUE,WAIC=TRUE,
  monitors=c('beta0','beta1','serr','sdsub','u'))
ranint.out$summary$all.chains
ranint.out$WAIC
library(coda)
gelman.diag(ranint.out$samples)
allsamps <- as.matrix(ranint.out$samples)
raftery.diag(allsamps)
effectiveSize(allsamps)
write.table(allsamps,file='ranintBayes.txt')

####################################
### Random Intercepts and Slopes ###
####################################

ranintslpCode <- nimbleCode({
  for (i in 1:N){
    response[i] ~ dnorm(mu[i],sd=serr)
    mu[i] <- beta0[tmtnum[i]] + beta1[tmtnum[i]]*ztime[i] + u0[subnum[i]] + u1[subnum[i]]*ztime[i]
  }
  for (i in 1:3){
    beta0[i] ~ dnorm(80,sd=5)
    beta1[i] ~ dnorm(0,sd=2)
  }
  for (i in 1:57){
    u0[i] ~ dnorm(0,sd=sdsub)
    u1[i] ~ dnorm(0,sd=sdslp)
  }
  serr ~ dgamma(1.1,1)
  sdsub ~ dgamma(1.1,.5)
  sdslp ~ dgamma(1.1,3)
})
N <- 399
tmtnum <- weight$tmtnum
subnum <- weight$subnum
response <- weight$response
ztime <- weight$ztime
ranintslpConsts <- list(N=N,tmtnum=tmtnum,subnum=subnum,ztime=ztime)
ranintslpData <- list(response=response)
ranintslp.out <- nimbleMCMC(
  code=ranintslpCode,
  constants = ranintslpConsts,
  data=ranintslpData,
  nchains=12,niter=102000,
  nburnin=2000,thin=100,
  samplesAsCodaMCMC = TRUE,
  summary=TRUE,WAIC=TRUE,
  monitors=c('beta0','beta1','serr','sdsub','sdslp','u0','u1'))
ranintslp.out$summary$all.chains
ranintslp.out$WAIC
gelman.diag(ranintslp.out$samples)
allsampsrc <- as.matrix(ranintslp.out$samples)
raftery.diag(allsampsrc)
effectiveSize(allsampsrc)
write.table(allsampsrc,file='ranintslpBayes.txt')

#############################
### Multivariate with AR1 ###
#############################

#first data must be reshaped so response is a vector
# of seven entries for each subject
weightord <- weight[order(weight$subnum),]  # order by subject
weightord <- weightord[order(weightord$tmtnum),] # order by tmt number
# put in matrix shape
dim(weightord)
#  number of subects is 399/7
resp <- t(matrix(weightord$response,7,57)) # row is response for each ind.
treatment <- t(matrix(weightord$tmtnum,7,57))  
weeknum <- t(matrix(weightord$ztime,7,57))
# First 20 rows are control 
# next 16 rows are RI
# last 21 rows are WI
# ar1 model
ar1Code <- nimbleCode({
  for (i in 1:nSub){
    y[i,1:7] ~ dmnorm(mu[i,1:7],cov=tau[1:7,1:7])
    for (j in 1:nRep){
      mu[i,j] <- b0[trt[i,j]]+b1[trt[i,j]]*tz[i,j]
    }
  }
  for (i in 1:nTmt){
    b0[i] ~ dnorm(80,sd=5)
    b1[i] ~ dnorm(0,sd=2)
  }
  rho ~ dbeta(1,1)
  s2err ~ dunif(0,50)
  for (i in 1:nRep){
    for (j in 1:nRep){
      tau[i,j] <- s2err*rho^(abs(i-j))
    }
  }
})
ar1Consts <- list(nSub=57,nRep=7,nTmt=3,
                  trt=treatment,
                  tz=weeknum)
ar1Data <- list(y=resp)
ar1.out <- nimbleMCMC(code=ar1Code,
                      data=ar1Data,
                      constants=ar1Consts,
                      nchains=5,niter=102000,
                      nburnin=2000,thin=50,
                      samplesAsCodaMCMC = TRUE,
                      summary=TRUE,WAIC=TRUE,
                      monitors=c('b0','b1','s2err','rho'))
ar1.out$summary$all.chains
ar1.out$WAIC
gelman.diag(ar1.out$samples)
allsampsar1 <- as.matrix(ar1.out$samples)
raftery.diag(allsampsar1)
effectiveSize(allsampsar1)
# Check the traceplot of the s2err because we used a uniform prior
# Make sure you fully span the space, not getting a 'haircut'
plot(allsampsar1[,8], type = 'l')

#######################################
### Fully Unstructured Multivariate ###
#######################################

mnormCode <- nimbleCode({
  for (i in 1:nSub){
    y[i,1:7] ~ dmnorm(mu[i,1:7],cov=tau[1:7,1:7])
    for (j in 1:nRep){
      mu[i,j] <- b0[trt[i,j]]+b1[trt[i,j]]*tz[i,j]
    }
  }
  for (i in 1:nTmt){
    b0[i] ~ dnorm(80, sd=5)
    b1[i] ~ dnorm(0,sd=2)
  }
  tau[1:7,1:7] ~ dinvwish(s[1:7,1:7],7)
  
})
s=matrix(6,7,7)
diag(s) <- 12
mnormConsts <- list(nSub=nSub,nRep=nRep,nTmt=nTmt,
                    trt=treatment,
                    tz=weeknum,
                    s=s)
mnormData <- list(y=resp)
mnorm.out <- nimbleMCMC(code=mnormCode,
                        data=mnormData,
                        constants=mnormConsts,
                        nchains=5,niter=12000,
                        nburnin=2000,thin=5,
                        samplesAsCodaMCMC = TRUE,
                        summary=TRUE,WAIC=TRUE,
                        monitors=c('b0','b1','tau'))
mnorm.out$summary$all.chains
mnorm.out$WAIC
gelman.diag(mnorm.out$samples[,1:6])
allsampsmv <- as.matrix(mnorm.out$samples)
raftery.diag(allsampsmv)
effectiveSize(allsampsmv)

#######################
######## TESTS ########
#######################

allsampsrc[1,]
# b01 = b03
# Are the intercepts between treatments the same? If they look different, we'd need to talk to 
# the designer to make sure it was designed/randomized correctly. If there is no difference, we probably
# wouldn't need to talk to them
plot(density(allsampsrc[,1]-allsampsrc[,3]))
mean(allsampsrc[,1]<allsampsrc[,3])

# Differences in slopes. Does one of the slopes appear significantly better?
plot(density(allsampsrc[,4]-allsampsrc[,6]))
mean(allsampsrc[,4]<allsampsrc[,6])
mean(allsampsrc[,5]<allsampsrc[,6])
plot(density(allsampsrc[,6]-allsampsrc[,5]))

#plot confidence intervals
# focus on WI
dim(allsampsrc) # 12000 rows
ztest <- 0:6
allsampsrc[1:12,6] # draws for the slope
test <- outer(allsampsrc[,6],ztest) # multiply all samples by all the zs
dim(test) # 12000 rows by 7 columns of time
ltest <- allsampsrc[,3] + test # add on the intercept
# Plot all 12k lines
plot(ztest,ltest[1,],type='l',col='lightgray',ylim=c(78,88))
for (i in 2:12000){lines(ztest,ltest[i,],type='l',col='lightgray')}
# Add the credible interval
lines(ztest,apply(ltest,2,mean),col='black')
lines(ztest,apply(ltest,2,quantile,.025),col='blue')
lines(ztest,apply(ltest,2,quantile,.975),col='blue')

# CI at different times
quantile(ltest[,1],c(.025,.975))
quantile(ltest[,2],c(.025,.975))
quantile(ltest[,3],c(.025,.975))
quantile(ltest[,4],c(.025,.975))
quantile(ltest[,5],c(.025,.975))
quantile(ltest[,6],c(.025,.975))
quantile(ltest[,7],c(.025,.975))

# prediction in treatment WI
# 7 is SD_slp, 8 is SD_sub, 9 is SD_err
p1test <- outer((allsampsrc[,6]+rnorm(12000,0,allsampsrc[,7])),ztest)
# Basically just drawing random errors
ptest <- (allsampsrc[,3]+rnorm(12000,0,allsampsrc[,8])) + rnorm(12000,0,allsampsrc[,9]) + p1test
dim(ptest)
plot(ztest,ptest[1,],type='l',col='lightgray',ylim=c(67,100))
for (i in 2:12000){lines(ztest,ptest[i,],type='l',col='lightgray')}
lines(ztest,apply(ptest,2,mean),col='black')
lines(ztest,apply(ptest,2,quantile,.025),col='blue')
lines(ztest,apply(ptest,2,quantile,.975),col='blue')

# PI at different times
quantile(ptest[,1],c(.025,.975))
quantile(ptest[,2],c(.025,.975))
quantile(ptest[,3],c(.025,.975))
quantile(ptest[,4],c(.025,.975))
quantile(ptest[,5],c(.025,.975))
quantile(ptest[,6],c(.025,.975))
quantile(ptest[,7],c(.025,.975))
