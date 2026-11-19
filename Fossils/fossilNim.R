fossil <- read.table('fossil.dat',header=T)
fossilS <- fossil[order(fossil$age),]
y <- fossilS$strontium.ratio*1000
library(nlme)
grp <- rep(1,106)  
fossilSG <- groupedData(strontium.ratio~age|grp,fossilS)
library(nimble)
knots <- default.knots(fossilSG$age,25)
knots
Z <- outer(fossilSG$age,knots,'-')
Z <- Z*(Z>0)
dim(Z)
Z[1:10,1]
Z[1:10,25]
Z[96:106,25]
smoothFossilCode <- nimbleCode({
  for (i in 1:106){
    y[i] ~ dnorm(mu[i],var=s2error)
    mu[i] <- b0 + b1*age[i] + u1*Z[i,1] + u2*Z[i,2] +
      u3*Z[i,3] + u4*Z[i,4] + u5*Z[i,5] +
      u6*Z[i,6] + u7*Z[i,7] + u8*Z[i,8] +
      u9*Z[i,9] + u10*Z[i,10] + u11*Z[i,11] +
      u12*Z[i,12] + u13*Z[i,13] + u14*Z[i,14] +
      u15*Z[i,15] + u16*Z[i,16] + u17*Z[i,17] +
      u18*Z[i,18] + u19*Z[i,19] + u20*Z[i,20] +
      u21*Z[i,21] + u22*Z[i,22] + u23*Z[i,23] +
      u24*Z[i,24] + u25*Z[i,25]
  }
  b0 ~ dnorm(704,sd=.5)
  b1 ~ dnorm(.04,sd=.01)
  u1 ~ dnorm(0,var=s2u)
  u2 ~ dnorm(0,var=s2u)
  u3 ~ dnorm(0,var=s2u)
  u4 ~ dnorm(0,var=s2u)
  u5 ~ dnorm(0,var=s2u)
  u6 ~ dnorm(0,var=s2u)
  u7 ~ dnorm(0,var=s2u)
  u8 ~ dnorm(0,var=s2u)
  u9 ~ dnorm(0,var=s2u)
  u10 ~ dnorm(0,var=s2u)
  u11 ~ dnorm(0,var=s2u)
  u12 ~ dnorm(0,var=s2u)
  u13 ~ dnorm(0,var=s2u)
  u14 ~ dnorm(0,var=s2u)
  u15 ~ dnorm(0,var=s2u)
  u16 ~ dnorm(0,var=s2u)
  u17 ~ dnorm(0,var=s2u)
  u18 ~ dnorm(0,var=s2u)
  u19 ~ dnorm(0,var=s2u)
  u20 ~ dnorm(0,var=s2u)
  u21 ~ dnorm(0,var=s2u)
  u22 ~ dnorm(0,var=s2u)
  u23 ~ dnorm(0,var=s2u)
  u24 ~ dnorm(0,var=s2u)
  u25 ~ dnorm(0,var=s2u)
  s2error ~ dgamma(2,500)
  s2u ~ dgamma(2,500)
})
age=fossilSG$age
smfData <- list(y=y)
smfConsts <- list(Z=Z,age=age)
smfinits <- list(b0=704,b1=0,s2error=.0006,s2u=.0004,u1=-.002,u2=-.006,
                 u3=-.02,u4=-.014,u5=.016,
                 u6=-.006,u7=-.009,u8=-.01,u9=-.003,
                 u10=0-.005,u11=-.008,u12=-.007,
                 u13=-.006,u14=.002,u15=.006,
                 u16=.015,u17=.006,u18=0.009,
                 u19=0.022,u20=0.013,u21=0.025,
                 u22=-.007,u23=-.026,u24=-.015,
                 u25=-.006)
smf.out <- nimbleMCMC(
  code = smoothFossilCode,
  data = smfData,
  constants = smfConsts,
  inits = smfinits,
  nchains = 5, niter = 2010000,
  nburnin=10000, thin=1000,
  samplesAsCodaMCMC = TRUE,
  WAIC=TRUE, summary = TRUE,
  monitors = c('b0','b1','s2error','s2u','u1','u2','u3','u4','u5',
               'u6','u7','u8','u9','u10','u11','u12','u13',
               'u14','u15','u16','u17','u18','u19','u20',
               'u21','u22','u23','u24','u25')
)
smf.out$summary
smf.out$WAIC
library(coda)
gelman.diag(smf.out$samples)
allsampsmf <- as.matrix(smf.out$samples)
raftery.diag(allsampsmf)
effectiveSize(allsampsmf)
X <- cbind(1,age)
beta <- allsampsmf[,1:2]
uhat <- allsampsmf[,c(5,16,22:28,6:15,17:22)]
fhatmn <- X%*%apply(beta,2,mean) + Z%*%apply(uhat,2,mean)
plot(fossilSG$age,y,ylim=c(707.2,707.5))
lines(age,fhatmn,col='blue')
