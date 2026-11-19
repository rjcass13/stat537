weight <- read.table('weight.dat',header=TRUE)
sub <- paste(substr(weight$tmt,1,1),sprintf("%02d",weight$subject),sep='')
weight$sub <- sub
head(weight)
ztime <- weight$time - 1
weight$ztime <- ztime
weightsort <- weight[order(weight$sub),]
weightsort[1:14,]
library(nlme)
WeightG <- groupedData(response ~ time | sub,data=weight)
plot(WeightG)
plot(WeightG,outer=~tmt)
fitlm <- lm(response~-1+as.factor(tmt)+ztime:as.factor(tmt),data=weight)
summary(fitlm)
logLik(fitlm)
library(lme4)
fitcs <- lmer(response ~-1+as.factor(tmt)+ztime:as.factor(tmt) + (1|sub),data=weight)
summary(fitcs)
logLik(fitcs)
fitcsg <- gls(response ~-1+as.factor(tmt)+ztime:as.factor(tmt),correlation = corCompSymm(form=~1|sub),data=weight)
logLik(fitcsg)
fitrc <- lmer(response ~ -1 + as.factor(tmt)+ztime:as.factor(tmt) + (ztime|sub),data=weight)
summary(fitrc)
logLik(fitrc)
fitrc1 <- lmer(response ~ -1 + as.factor(tmt)+ztime:as.factor(tmt) + (ztime||sub),data=weightsort)
summary(fitrc1)
logLik(fitrc1)
library(lmerTest)
fitTest <- lmerTest::lmer(response ~ as.factor(tmt)+ztime:as.factor(tmt) + (ztime||sub),data=weightsort)
anova(fitTest)
bhat <- t(t(fixef(fitrc1)))
bhat
uhat <- ranef(fitrc1)$sub
uhat
test <- matrix(as.matrix(uhat),114,1)
uhat <- test
uhat
X <- model.matrix(~-1+as.factor(tmt)+ztime:as.factor(tmt),data=weightsort)
dim(X)
Z <- model.matrix(response~-1+as.factor(sub)+ztime:sub,data=weightsort)
dim(Z)
Z[2,]
summary(fitrc1)
R <- .6337*diag(1,399)
g <- rbind(c(8.8853,0),c(0,.1234))
g
G <- g%x%diag(1,57)
dim(G)
G[1:7,1:7]
G[56:60,56:60]
sigma <- Z%*%G%*%t(Z)+R
sigma[1:7,1:7]
v11 <- solve(t(X)%*%solve(sigma)%*%X)
v22 <- solve(t(Z)%*%R%*%Z + solve(G)) + G%*%t(Z)%*%solve(sigma)%*%X%*%v11%*%t(X)%*%solve(sigma)%*%Z%*%G
v12 <- -v11%*%t(X)%*%solve(sigma)%*%Z%*%G
v21 <- t(v12)
chat <- rbind(cbind(v11,v12),cbind(v21,v22))
z56 <- rep(0,56)
z56 <- t(z56)
csub1 <- cbind(1,0,0,6,0,0,1,z56,6,z56)
csub1
bhatuhat <- rbind(bhat,uhat)
bhatuhat
csub1%*%bhatuhat
# confidence interval for subject 1 in control at ztime 6
# actual response is 87
csub1%*%bhatuhat + qt(.975,55)*sqrt(csub1%*%chat%*%t(csub1))
csub1%*%bhatuhat - qt(.975,55)*sqrt(csub1%*%chat%*%t(csub1))
#prediction interval (which is a little silly here)
csub1%*%bhatuhat + qt(.975,55)*sqrt(csub1%*%chat%*%t(csub1)+.6337)
csub1%*%bhatuhat - qt(.975,55)*sqrt(csub1%*%chat%*%t(csub1)+.6337)
vhatbhat <- chat[1:6,1:6]
cw6 <- rbind(c(0,0,1,0,0,6))
cw6
ucl <- cw6%*%bhat + qt(.975,55)*sqrt(cw6%*%vhatbhat%*%t(cw6))
lcl <- cw6%*%bhat - qt(.975,55)*sqrt(cw6%*%vhatbhat%*%t(cw6))
c(lcl,cw6%*%bhat,ucl)
cw <- cbind(0,0,1,0,0,0:6)
cw
plot(1:7,cw%*%bhat,xlim=c(0,8),ylim=c(75,90),type='l')
ucl <- cw%*%bhat + qt(.975,55)*sqrt(diag(cw%*%vhatbhat%*%t(cw)))
lcl <- cw%*%bhat - qt(.975,55)*sqrt(diag(cw%*%vhatbhat%*%t(cw)))
cbind(lcl,cw%*%bhat,ucl)
lines(1:7,ucl,col='green')
lines(1:7,lcl,col='green')
cr <- cbind(0,1,0,0,0:6,0)
lines(1:7,cr%*%bhat,lty=3,col='blue')
ucl <- cr%*%bhat + qt(.975,55)*sqrt(diag(cr%*%vhatbhat%*%t(cr)))
lcl <- cr%*%bhat - qt(.975,55)*sqrt(diag(cr%*%vhatbhat%*%t(cr)))
cbind(lcl,cr%*%bhat,ucl)
lines(1:7,ucl,col='orange',lty=3)
lines(1:7,lcl,col='orange',lty=3)
crw1 <- cbind(0,-1,1,0,0,0)
tstat <- (crw1%*%bhat)/sqrt(crw1%*%vhatbhat%*%t(crw1))
tstat
2*(1-pt(tstat,55))
crw6 <- cbind(0,-1,1,0,-6,6)
tstat <- (crw6%*%bhat)/sqrt(crw6%*%vhatbhat%*%t(crw6))
tstat
2*(1-pt(tstat,55))
csrw <- cbind(0,0,0,0,-1,1)
tstat <- (csrw%*%bhat)/sqrt(csrw%*%vhatbhat%*%t(csrw))
tstat
2*(1-pt(tstat,55))
# This next section is illustrative, but incorrect
upl <- cr%*%bhat + qt(.975,55)*sqrt(diag(cr%*%vhatbhat%*%t(cr)+.6337))
lpl <- cr%*%bhat - qt(.975,55)*sqrt(diag(cr%*%vhatbhat%*%t(cr)+.6337))
lines(1:7,upl,col='blue',lty=3)
lines(1:7,lpl,col='blue',lty=3)
upl <- cw%*%bhat + qt(.975,55)*sqrt(diag(cw%*%vhatbhat%*%t(cw)+.6337))
lpl <- cw%*%bhat - qt(.975,55)*sqrt(diag(cw%*%vhatbhat%*%t(cw)+.6337))
lines(1:7,upl,col='blue',lty=1)
lines(1:7,lpl,col='blue',lty=1)
plot(1:7,cr%*%bhat,xlim=c(0,8),ylim=c(70,90),type='l',col='blue')
# This is the correct way to do prediction intervals
z <- cbind(1,0:6)
r <- .6337*diag(1,7)
vardata <- z%*%g%*%t(z)+r
upl <- cr%*%bhat + qt(.975,55)*sqrt(diag(cr%*%vhatbhat%*%t(cr)+vardata))
lpl <- cr%*%bhat - qt(.975,55)*sqrt(diag(cr%*%vhatbhat%*%t(cr)+vardata))
lines(1:7,upl,col='blue',lty=3)
lines(1:7,lpl,col='blue',lty=3)
# other potential models
fitar1 <- gls(response ~-1+as.factor(tmt)+ztime:as.factor(tmt),correlation = corAR1(form=~1|sub),data=weight)
summary(fitar1)
3.282^2
logLik(fitar1)
getVarCov(fitar1)
vhatbhat <- vcov(fitar1)
vhatbhat
bb <- matrix(0,7,7)
bb
for (i in 1:7){for (j in 1:7){bb[i,j]<- ifelse(i==j,10.77,10.77*.95^(abs(i-j)))}}
bb
fitar1w <- gls(response ~-1+as.factor(tmt)+ztime:as.factor(tmt),correlation = corAR1(form=~1|sub),weights=varIdent(form=~1|tmt),data=weight)
summary(fitar1w)
logLik(fitar1w)
anova(fitar1,fitar1w)
fitar1wzt <- gls(response ~-1+as.factor(tmt)+ztime:as.factor(tmt),correlation = corAR1(form=~1|sub),weights=varIdent(form=~1|ztime),data=weight)
anova(fitar1,fitar1wzt)
summary(fitar1wzt)
logLik(fitar1wzt)
getVarCov(fitar1wzt)
fitun <- gls(response ~ -1+as.factor(tmt)+ztime:as.factor(tmt),correlation=corSymm(form=~1|sub),weights=varIdent(form=~1|ztime),data=weight)
logLik(fitun)
getVarCov(fitun)
anova(fitar1,fitun)
summary(fitun)
AIC(fitrc1)
AIC(fitar1)
AIC(fitun)
AIC(fitar1w)
AIC(fitar1wzt)
BIC(fitrc1)
BIC(fitar1)
BIC(fitun)
BIC(fitar1w)
BIC(fitar1wzt)
keep <- sample(1:399,300)
weightun <- weight[keep,]
dim(weightun)
fitrc1u <- lmer(response ~ -1 + as.factor(tmt)+ztime:as.factor(tmt) + (ztime||sub),data=weightun)
fitar1u <- gls(response ~-1+as.factor(tmt)+ztime:as.factor(tmt),correlation = corAR1(form=~1|sub),data=weightun)
fitar1wu <- gls(response ~-1+as.factor(tmt)+ztime:as.factor(tmt),correlation = corAR1(form=~1|sub),weights=varIdent(form=~1|tmt),data=weightun)
AIC(fitrc1u)
AIC(fitar1u)
AIC(fitar1wu)
BIC(fitrc1u)
BIC(fitar1u)
BIC(fitar1wu)
