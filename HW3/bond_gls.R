# bond example using gls
dat <- read.table('bond.dat')
names(dat) <- c('ingot','metal','strength')
dat
W <- model.matrix(~-1 + as.factor(metal),dat)
W
olsbhat <- solve(t(W)%*%W)%*%t(W)%*%dat$strength
olsbhat
sse <- t(dat$strength-W%*%olsbhat)%*%(dat$strength-W%*%olsbhat)
sse
dim(W)
mse <- as.numeric(sse)/18 # 18 = 21 - 3 = n-p
volsbhat <- mse*solve(t(W)%*%W)
volsbhat
cm <- rbind(c(1,-1,0),c(0,1,-1))
cm
cm%*%olsbhat
# The ftest here is basically looking at how much of the difference is due to 
# noise versus actual signal. 
fhat <- t(cm%*%olsbhat)%*%solve(cm%*%volsbhat%*%t(cm))%*%cm%*%olsbhat/2
fhat
1-pf(fhat,2,18)

# Now a different way
# GLS: Generalized Least Squares
library(nlme)
# corCompSymm: Compound Symmetric, this defines the covariance structure. 
# This makes each ingot have its own covariance, as described in the notes
fit1 <- gls(strength ~ -1 + as.factor(metal),correlation=corCompSymm(form=~1|ingot),method='ML',data=dat)
summary(fit1)
# The penalty is -2*number of effects + number of variance effects, in our case sigmaI and sigmaE
# You can get the AIC and BIC from the log likelihood: -2l + penalties
# Standardized residuals we expect to be more or less between -2 and 2

# the covariance structure
# total variance down diagonal
4.324617^2
# off diagonal covariance (ingot to ingot)
18.70231*.5246614 # (.5... is the rho value)
# pure error variance (sigma^2_e)
18.70231 - 18.70231*.5246614


# wrong test
anova(fit1) # original fit is mu1 = mu2 = mu3 = 0 because we removed the intercept

#correct test (maybe)
fit2 <- gls(strength ~ as.factor(metal),correlation=corCompSymm(form=~1|ingot),method='ML',data=dat)
anova(fit2)
fit2 <- gls(strength ~ as.factor(metal),correlation=corCompSymm(form=~1|ingot),method='REML',data=dat)
anova(fit2)

#do F by first principles - safer
sigma <- matrix(.5246614,3,3)
sigma
diag(sigma) <- c(1,1,1)
sigma
sigma <- 4.324617^2*sigma
sigma
sigma <- kronecker(diag(7),sigma) # Similar to 'outer' joins, makes a 7x7 identity, then populates 
# each diagonal with sigma
sigma
sigmai <- solve(sigma)

W <- model.matrix(~-1 + as.factor(metal),dat)
bhat <- solve(t(W)%*%sigmai%*%W)%*%t(W)%*%sigmai%*%dat$strength
# Compare the calculated estimates
fit1$coefficients
bhat

varbhat <- solve(t(W)%*%sigmai%*%W)
varbhat

cm <- rbind(c(1,-1,0),c(0,1,-1))
cm
cm%*%bhat
fhat <- t(cm%*%bhat)%*%solve(cm%*%varbhat%*%t(cm))%*%cm%*%bhat/2
fhat
# MLEs of variances are generaly biased too low
# This means the F is generally too big
# Even if it's unbiased, it's still an approximate F because there is no closed form 
# solution for the variance. We're assuming these estimates for variances are correct,
# but they're not exactly the true value
1-pf(fhat,2,18)
1-pf(fhat,2,16)

# Likelihood Ratio Test
# Test that there's a difference in metal strength
fit3 <- gls(strength ~ 1,correlation=corCompSymm(form=~1|ingot),method='ML',data=dat)
anova(fit1,fit3)
# Calculate the p-value directly
pchisq(10.11648, 2, lower.tail = FALSE)

# Let's try a general covariance structure
fit4 <- gls(strength ~ -1+as.factor(metal),correlation=corSymm(form=~1|ingot),weights=varIdent(form=~1|metal),method='ML',data=dat)
summary(fit4)

# Variances as pulled from Summary: the diagonals of the covariance matrix
3.940268^2 # Metal 1 Variance (Residual Standard Error squared)
(1.4421742*3.940268)^2 # Metal 2 variance (RSE * ratio provided in 'variance function' section)
(.7307121*3.940268)^2 # Metal 3 Variance (RSE * ratio provided in 'variance function' section)

# Covariances as pulled from summary: the off diagonals for the covariance matrix
# rho_{12} = sigma_{12} / (sigma_1 * sigma_2) => sigma_{12} = rho_{12}*sigma_1*sigma_2
0.909*(3.940268)*((1.4421742*3.940268))
0.384*(3.940268)*((.7307121*3.940268))
0.290*(1.4421742*3.940268)*((.7307121*3.940268))

# Full covariance matrix
sigma <- rbind(c(15.52,20.35,4.36),c(20.35,32.29,4.74),c(4.36,4.74,8.29))
sigma
sigma <- kronecker(diag(7),sigma)
sigma
sigmai <- solve(sigma)
bhat <- solve(t(W)%*%sigmai%*%W)%*%t(W)%*%sigmai%*%dat$strength
bhat
varbhat <- solve(t(W)%*%sigmai%*%W)
varbhat
cm <- rbind(c(1,-1,0),c(0,1,-1))
cm
cm%*%bhat
fhat <- t(cm%*%bhat)%*%solve(cm%*%varbhat%*%t(cm))%*%cm%*%bhat/2
fhat
1-pf(fhat,2,18)
1-pf(fhat,2,9)

# test metal significance using LRT
fit5 <- gls(strength ~ 1,correlation=corSymm(form=~1|ingot),weights=varIdent(form=~1|metal),method='ML',data=dat)
summary(fit5)
anova(fit4,fit5)

# test covariance structure using LRT
anova(fit4,fit1)
# With 4 (9-5) degrees of freedom, fit4 increased the likelihood (less negative) as compared to fit1

# covariance structure using AIC
c(AIC(fit1),AIC(fit4))
# Want AIC to be low, so still indicates fit4

# covariance structure using BIC
c(BIC(fit1),BIC(fit4))
# Want BIC to be low, so still indicates fit4