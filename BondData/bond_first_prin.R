bond <- read.table('bond.dat')
names(bond) <- c('ingot','metal','strength')

z <- model.matrix(~-1+as.factor(bond$ingot))
z
diag(1,3)
rep(1,7)
t(t(rep(1,7)))
x <- kronecker(t(t(rep(1,7))),diag(1,3))
x
r <-diag(10,21)
r
g <- diag(12,7)
g
z%*%g%*%t(z)
sigma <- z%*%g%*%t(z) + r
sigma
y <- t(t(bond$strength))
y

# These estimates are unbiased, but do not correctly account for the variance
bhat1 <- solve(t(x)%*%x)%*%t(x)%*%y
bhat1
si <- solve(sigma)
si
bhat2 <- solve(t(x)%*%si%*%x)%*%t(x)%*%si%*%y
bhat2

n <- length(y)
n
ss1 <- t(y-x%*%bhat1)%*%(y-x%*%bhat1)
ss1
p <- dim(x)[2]
p
ms <- ss1/(n-p)
ms
as.numeric(ms)
ms <- as.numeric(ms)

# This is the correct error
s2bhat1 <- ms*solve(t(x)%*%x)
s2bhat1
s2bhat2 <- solve(t(x)%*%si%*%x)
s2bhat2

# First Row: Mean Metal1 = Mean Metal2
# Second Row: Mean Metal2 = Mean Metal3
cm <- rbind(c(1,-1,0),c(1,0,-1))
cm

# What's the degrees of freedom?
# There's not really a 'correct' one, but use one that you can justify
fstat1 <- t(cm%*%bhat1)%*%solve(cm%*%s2bhat1%*%t(cm))%*%cm%*%bhat1/2
fcrit1 <- qf(.95,2,18)
fstat1
fcrit1
1-pf(as.numeric(fstat1),2,18)
fstat2 <- t(cm%*%bhat2)%*%solve(cm%*%s2bhat2%*%t(cm))%*%cm%*%bhat2/2
fcrit2 <- qf(.95,2,14)
fstat2
fcrit2
1-pf(as.numeric(fstat2),2,14)

#generate random data with the covariance structure in sigma
sigma
s <- chol(sigma)
t(s)%*%s
egen <- rnorm(21)
egen
enew <- t(s)%*%egen
enew
# now add the mean structure
ynew <-rep(c(70,75,70),7)+enew
ynew
solve(t(x)%*%si%*%x)%*%t(x)%*%si%*%ynew

# Using the MASS library to generate test data
library(MASS)
mvrnorm(1,rep(c(70,75,70),7),sigma)
# Compare to the previous ynews
cbind(ynew,mvrnorm(1,rep(c(70,75,70),7),sigma))

# Another method to generate data
u <- rnorm(7,0,sqrt(12))
u
e <- rnorm(21,0,sqrt(10))
e
# first obs would be
70 + u[1] +e[1]
# second obs
75 + u[1] + e[2]
# fourth obs
70 + u[2] + e[4]
