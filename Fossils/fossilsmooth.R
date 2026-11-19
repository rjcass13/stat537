
fossil <- read.table('fossil.dat',header=T)
fossilS <- fossil[order(fossil$age),]
age <- fossilS$age
agec <- (age-mean(age))/sd(age)
fossilS$agec <- agec
y <- fossilS$strontium.ratio*1000
plot(fossil$age,fossil$strontium.ratio)
library(nlme)
grp <- rep(1,106)  
fossilSG <- groupedData(strontium.ratio~age|grp,fossilS)
#knots <- default.knots(fossilSG$age,25)
#knots <- seq(94,121,length=25)
#knots <- seq(93,122,length=200)
knots <- default.knots(fossilSG$age,100)
knots
Z <- outer(fossilSG$age,knots,'-')
Z <- Z*(Z>0)
dim(Z)
Z[1:10,1]
Z[15:20,1]
Z[,1]
Z[,2]
fit.ex <- lm(strontium.ratio~age +Z,data=fossilSG)
yhat <- fitted(fit.ex)
plot(fossilSG$age,fossilSG$strontium.ratio)
lines(fossilSG$age,yhat)
fit.fossil <- lme(strontium.ratio~age,data=fossilSG,random=pdIdent(~-1+Z))
beta <- fit.fossil$coef$fixed
beta
uhat <- unlist(fit.fossil$coef$random)
uhat
plot(fossilSG$age,fossilSG$strontium.ratio,ylim=c(.70715,.70755))
X <- cbind(1,fossilSG$age)
fhat <- X%*%beta + Z%*%uhat
#lines(fossilSG$age,yhat,col='red')
lines(fossilSG$age,fhat,col='blue',lwd=2)
sigma.e <- fit.fossil$sigma
sigma.u <- sigma.e*exp(unlist(fit.fossil$modelStruct))
s2e <- sigma.e^2
s2e
s2u <- sigma.u^2
s2u
r <- diag(s2e,dim(Z)[1])
g <- diag(s2u,dim(Z)[2])
v <- Z%*%g%*%t(Z) + r
vcov(fit.fossil)
solve(t(X)%*%solve(v)%*%X)
ri <- solve(r)
gi <- solve(g)
c11 <- t(X)%*%ri%*%X
c12 <- t(X)%*%ri%*%Z
c21 <- t(Z)%*%ri%*%X
c22 <- t(Z)%*%ri%*%Z + gi
ci <- rbind(cbind(c11,c12),cbind(c21,c22))
chat <- solve(ci)
dim(chat)
XB <- cbind(X,Z)
uci <- X%*%beta+Z%*%uhat + qt(1-.05/2,dim(X)[1]-dim(X)[2])*sqrt(diag(XB%*%chat%*%t(XB)))
lines(fossilSG$age,uci,lty=2,col='red')
lci <- X%*%beta+Z%*%uhat - qt(1-.05/2,dim(X)[1]-dim(X)[2])*sqrt(diag(XB%*%chat%*%t(XB)))
lines(fossilSG$age,lci,lty=2,col='red')
upi <- X%*%beta+Z%*%uhat + qt(1-.05/2,dim(X)[1]-dim(X)[2])*sqrt(diag(XB%*%chat%*%t(XB)+s2e))
lines(fossilSG$age,upi,lty=3,col='green')
lpi <- X%*%beta+Z%*%uhat - qt(1-.05/2,dim(X)[1]-dim(X)[2])*sqrt(diag(XB%*%chat%*%t(XB)+s2e))
lines(fossilSG$age,lpi,lty=3,col='green')
legend('bottomleft',c('Smooth','Conf. Int.','Pred. Int.'),col=c('blue','red','green'),lty=c(1,2,3))


# example of cubic splines
# Show the difference between 20 and 50 knots with this example.


x2 <- fossilSG$age^2
x3 <- fossilSG$age^3
X <- cbind(1,fossilSG$age,x2,x3)
knots <- default.knots(fossilSG$age,25)
 #knots <- default.knots(fossilSG$age,50)
# knots <- default.knots(fossilSG$age,5)
# knots <- seq(94,121,length=50)
# knots <- seq(95,120,5)
Z <- outer(fossilSG$age,knots,'-')
Z <- Z*(Z>0)
Z <- Z^3
fit.fossil.cubic <- lme(strontium.ratio~-1+X,random=pdIdent(~-1+Z),fossilSG)
betac <- t(t(fit.fossil.cubic$coef$fixed))
uhatc <- t(t(unlist(fit.fossil.cubic$coef$random)))
fhatc <- X%*%betac + Z%*%uhatc
plot(fossilSG$age,fossilSG$strontium.ratio,ylim=c(.7072,.7075))
lines(fossilSG$age,fhatc,col='blue',lwd=2)


# example of cubic radial basis functions pp 35 Ngo & Wand

X <- cbind(1,fossilSG$age)
knots <- default.knots(fossilSG$age,25)
#knots <- seq(93,121,length=50)
svd.Omega <- svd(abs(outer(knots,knots,'-'))^3)
matrix.sqrt.Omega <- t(svd.Omega$v%*%(t(svd.Omega$u)*sqrt(svd.Omega$d)))
Z <- t(solve(matrix.sqrt.Omega,t(abs(outer(fossilSG$age,knots,'-')^3))))
fit.fossil <- lme(strontium.ratio~-1+X,random=pdIdent(~-1+Z),fossilSG)
betar <- t(t(fit.fossil.radial$coef$fixed))
uhatr <- t(t(unlist(fit.fossil.radial$coef$random)))
fhatr <- X%*%betar + Z%*%uhatr
plot(fossilSG$age,fossilSG$strontium.ratio,ylim=c(.7072,.7075))
lines(fossilSG$age,fhatr,col='green',lwd=2)

# example of cubic splines with age centered

fossilSGagec <- groupedData(strontium.ratio~agec|grp,fossilS)
x2 <- agec^2
x3 <- agec^3
X <- cbind(1,agec,x2,x3)
#knots <- default.knots(fossilSG$age,10)
#knots <- default.knots(fossilSG$age,25)
knots <- default.knots(fossilSGagec$agec,50)
# knots <- seq(94,121,length=50)
# knots <- seq(95,120,5)
Z <- outer(fossilSGagec$agec,knots,'-')
Z <- Z*(Z>0)
Z <- Z^3
fit.fossil.cubic <- lme(strontium.ratio~-1+X,random=pdIdent(~-1+Z),fossilSGagec)
betac <- t(t(fit.fossil.cubic$coef$fixed))
uhatc <- t(t(unlist(fit.fossil.cubic$coef$random)))
fhatc <- X%*%betac + Z%*%uhatc
plot(fossilSGagec$agec,fossilSGagec$strontium.ratio,ylim=c(.7072,.7075))
lines(fossilSGagec$agec,fhatc,col='blue',lwd=2)
sigma.e <- fit.fossil.cubic$sigma
sigma.u <- sigma.e*exp(unlist(fit.fossil.cubic$modelStruct))
s2e <- sigma.e^2
s2e
s2u <- sigma.u^2
s2u
r <- diag(s2e,dim(Z)[1])
g <- diag(s2u,dim(Z)[2])
v <- Z%*%g%*%t(Z) + r
vcov(fit.fossil.cubic)
solve(t(X)%*%solve(v)%*%X)
ri <- solve(r)
gi <- solve(g)
c11 <- t(X)%*%ri%*%X
c12 <- t(X)%*%ri%*%Z
c21 <- t(Z)%*%ri%*%X
c22 <- t(Z)%*%ri%*%Z + gi
ci <- rbind(cbind(c11,c12),cbind(c21,c22))
chat <- solve(ci)
XB <- cbind(X,Z)
uci <- X%*%betac+Z%*%uhatc + qt(1-.05/2,dim(X)[1]-dim(X)[2])*sqrt(diag(XB%*%chat%*%t(XB)))
lines(fossilSGagec$agec,uci,lty=2,col='red')
lci <- X%*%betac+Z%*%uhatc - qt(1-.05/2,dim(X)[1]-dim(X)[2])*sqrt(diag(XB%*%chat%*%t(XB)))
lines(fossilSGagec$agec,lci,lty=2,col='red')
upi <- X%*%betac+Z%*%uhatc + qt(1-.05/2,dim(X)[1]-dim(X)[2])*sqrt(diag(XB%*%chat%*%t(XB)+s2e))
lines(fossilSGagec$agec,upi,lty=3,col='green')
lpi <- X%*%betac+Z%*%uhatc - qt(1-.05/2,dim(X)[1]-dim(X)[2])*sqrt(diag(XB%*%chat%*%t(XB)+s2e))
lines(fossilSGagec$agec,lpi,lty=3,col='green')

# example of cubic splines with orthonormal polynomials
source('orpol.R')
X <- orpol(agec,3)
dim(X)
head(X)
knots <- default.knots(X[,2],50)
# knots <- seq(94,121,length=50)
# knots <- seq(95,120,5)
Z <- outer(X[,2],knots,'-')
Z <- Z*(Z>0)
Z <- Z^3
fit.fossil.orth <- lme(strontium.ratio~-1+X,random=pdIdent(~-1+Z),fossilSG)
betao <- t(t(fit.fossil.orth$coef$fixed))
uhato <- t(t(unlist(fit.fossil.orth$coef$random)))
fhato <- X%*%betao + Z%*%uhato
plot(X[,2],fossilSGagec$strontium.ratio,ylim=c(.7072,.7075))
lines(X[,2],fhato,col='blue',lwd=2)
sigma.e <- fit.fossil.orth$sigma
sigma.u <- sigma.e*exp(unlist(fit.fossil.orth$modelStruct))
s2e <- sigma.e^2
s2e
s2u <- sigma.u^2
s2u
r <- diag(s2e,dim(Z)[1])
g <- diag(s2u,dim(Z)[2])
v <- Z%*%g%*%t(Z) + r
vcov(fit.fossil.orth)
solve(t(X)%*%solve(v)%*%X)
ri <- solve(r)
gi <- solve(g)
c11 <- t(X)%*%ri%*%X
c12 <- t(X)%*%ri%*%Z
c21 <- t(Z)%*%ri%*%X
c22 <- t(Z)%*%ri%*%Z + gi
ci <- rbind(cbind(c11,c12),cbind(c21,c22))
chat <- solve(ci)
XB <- cbind(X,Z)
uci <- X%*%betao+Z%*%uhato + qt(1-.05/2,50)*sqrt(diag(XB%*%chat%*%t(XB)))
lines(X[,2],uci,lty=2,col='red')
lci <- X%*%betao+Z%*%uhato - qt(1-.05/2,50)*sqrt(diag(XB%*%chat%*%t(XB)))
lines(X[,2],lci,lty=2,col='red')
upi <- X%*%betao+Z%*%uhato + qt(1-.05/2,50)*sqrt(diag(XB%*%chat%*%t(XB)+s2e))
lines(X[,2],upi,lty=3,col='green')
lpi <- X%*%betao+Z%*%uhato - qt(1-.05/2,50)*sqrt(diag(XB%*%chat%*%t(XB)+s2e))
lines(X[,2],lpi,lty=3,col='green')

### test for numerical sensitivity

fossilS <- fossil[order(fossil$age),]
y <- fossilS$strontium.ratio*1000
knots <- default.knots(fossilS$age,25)
knots
Z <- outer(fossilS$age,knots,'-')
Z <- Z*(Z>0)
fit.fossil.test <- lme(y~age,data=fossilSG,random=pdIdent(~-1+Z))
beta <- fit.fossil.test$coef$fixed
beta
uhat <- unlist(fit.fossil.test$coef$random)
uhat
plot(fossilSG$age,y,ylim=c(707.2,707.5))
X <- cbind(1,fossilS$age)
fhat <- X%*%beta + Z%*%uhat
lines(fossilSG$age,fhat,col='blue',lwd=2)
sigma.e <- fit.fossil.test$sigma
sigma.u <- sigma.e*exp(unlist(fit.fossil.test$modelStruct))
s2e <- sigma.e^2
s2e
s2u <- sigma.u^2
s2u
AIC(fit.fossil.test)
vcov(fit.fossil.test)
uhat
