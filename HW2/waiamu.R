#######################
#### In Class Code ####
#######################

adv <- read.table('advert.dat')
colnames(adv) <- c('salesvol', 'price', 'adv', 'tmt')
adv$tmtf <- as.factor(adv$tmt)
y <- adv$salesvol
# Removes the intercept
w <- model.matrix(y~tmtf-1,data=adv)
w
x <- model.matrix(y ~ tmtf,data=adv)
x

# Standard y representation: y = xB
# Can also be written as: Y = w mu
# We can also include a matrix and its inverse: y = w Ai A mu
# Where the following are true: w Ai = x, A mu = B

# So the 'A' matrix is the linear combinations of the mus we are working with

# This creates the A-inverse matrix using the x matrix, basically picking all the unique rows
ai <- rbind(
  cbind(1,0,0,0,0,0),
  cbind(1,1,0,0,0,0),
  cbind(1,0,1,0,0,0),
  cbind(1,0,0,1,0,0),
  cbind(1,0,0,0,1,0),
  cbind(1,0,0,0,0,1))

ai
w%*%ai
a <- solve(ai)
a
# Beta is A mu, so B1 = mu1, B2 = mu2 - mu1, B3 = mu3-mu1, etc.

# The actual Beta values: (x'x)^-1 x'y
solve(t(x)%*%x)%*%t(x)%*%y
# Now the mus using the w matrix
muhat <- solve(t(w)%*%w)%*%t(w)%*%y
# Again, Beta is A mu, so B1 = mu1, B2 = mu2 - mu1, B3 = mu3-mu1, etc.


# so the beta's are just linear combinations of the cell means
# the standard test of no difference in the means would be
anova(lm(y ~ tmtf,data=adv))
# how do we get the F of 8.4593?
# the last 5 beta's are what is needed
beta <- solve(t(x)%*%x)%*%t(x)%*%y
beta
sse <- t(y-x%*%beta)%*%(y-x%*%beta)
sse
sse <- as.numeric(sse)
mse <- sse/(18-6) # n-p (degrees of freedom), the p degrees of freedom 'used up' are the cell-means
mse
varbeta <- mse*solve(t(x)%*%x)
varbeta

# This is the F-test: (CB)'(CVar(B)C')CB/rank(C)
# In this case, since we got the x matrix, C is just the identity
t(beta[2:6])%*%solve(varbeta[2:6,2:6])%*%beta[2:6]/5
# so the f statistic using the beta's can be computed without
# worrying about the linear combinations since the formation of
# the x matrix took care of that
# Can we get the f statistic from the cell means model?
# Yes, with the right contrast matrix
cmat <- rbind(cbind(1,-1,0,0,0,0),cbind(1,0,-1,0,0,0),cbind(1,0,0,-1,0,0),cbind(1,0,0,0,-1,0),cbind(1,0,0,0,0,-1))
muhat
varmuhat <- mse*solve(t(w)%*%w)
t(cmat%*%muhat)%*%solve(cmat%*%varmuhat%*%t(cmat))%*%cmat%*%muhat/5
# so the advantage to the x matrix in computer programs is that
# they can use an identity for the linear contrasts of the 
# cell means, because the beta's are the appropriate linear contrasts

# now let's look at this as a 3 (price) x 2 (adv) design
# so the first cell mean is the 1,1 cell, the second cell mean
# is the 1,2 cell, the third cell mean is the 2,1 cell, etc.
#     | Radio | Newspaper |
# .59 |   1   |     2     |
# .60 |   3   |     4     |
# .64 |   5   |     6     |
# so the linear contrasts of the cell means that we want would be
# the overall mean, the contrast of price 1 with price 2, the contrast
# of price 2 with price 3, the contrast of radio with newspaper, and the
# two interaction degrees of freedom
a <- rbind(
  c(1/6,1/6,1/6,1/6,1/6,1/6),
  c(1/2,1/2,-1/2,-1/2,0,0),
  c(0,0,1/2,1/2,-1/2,-1/2),
  c(1/3,-1/3,1/3,-1/3,1/3,-1/3))
a
# we have one row for the mean, 2 rows for the effect of price, 1 row for
# the effect of advertising, and we need 2 rows for interaction.  We get
# these two rows as demonstrated in class
cint <- rbind(
  c(1,-1,-1,1,0,0),
  c(0,0,1,-1,-1,1))
a <- rbind(a,cint)
a
# first do the analysis the way you've been taught
anova(lm(salesvol ~ as.factor(price)*adv,data=adv))
# now with the X matrix we create as w%*%ai
ai <- solve(a)
X <- w%*%ai
X
betahat <- solve(t(X)%*%X)%*%t(X)%*%y
betahat
# test of price
sse <- t(y-X%*%betahat)%*%(y-X%*%betahat)
mse <- as.numeric(sse)/(18-6)
varbeta <- mse*solve(t(X)%*%X)
t(betahat[2:3])%*%solve(varbeta[2:3,2:3])%*%betahat[2:3]/2
# test of adv
betahat[4]%*%solve(varbeta[4,4])%*%betahat[4]/1
# test of interaction
t(betahat[5:6])%*%solve(varbeta[5:6,5:6])%*%betahat[5:6]/2
# now with the cell means
cp <- a[2:3,]
cp
t(cp%*%muhat)%*%solve(cp%*%varmuhat%*%t(cp))%*%cp%*%muhat/2
ca <- a[4,]
t(ca%*%muhat)%*%solve(ca%*%varmuhat%*%ca)%*%ca%*%muhat/1
ci <- a[5:6,]
t(ci%*%muhat)%*%solve(ci%*%varmuhat%*%t(ci))%*%ci%*%muhat/2


# now what if a cell is missing?
# let's remove all the data in the 1,1 cell
advmiss <- adv[c(-1,-2,-3),]
advmiss
# we can treat this as a one-way problem with 5 cells
w <- model.matrix(salesvol ~ -1 + tmtf,data=advmiss)
w
# the first column is all 0's since there is no 1,1 cell
# now solve(t(w)%*%w) won't work because of the first column
solve(t(w)%*%w)
# so let's remove the first column and work with the final 5 cols
w1 <- w[,2:6]
w1
ym <- advmiss$salesvol
muhat <- solve(t(w1)%*%w1)%*%t(w1)%*%ym
muhat
# now the trick is to get an X matrix that we can use
# since I took out the 1,1 cell we can use the same a and ai as before
# the w matrix and ai we get
Xm <- w%*%ai
Xm
# looking at Xm we see the first row of ai is missing
ai
# this is as it should be because there are no data in the 1,1 cell
# Now, while we have at Xm matrix, it is not full rank, so
# we can't use it to get estimates
solve(t(Xm)%*%Xm)
# It turns out the Xm matrix is of rank 5, so it would be possible to
# get a full rank matrix if we are thoughtful about the process.
# Since there are no data in the 1,1 cell, the first interaction df is
# not estimable.  IF we are willing to assume this interaction df is
# exactly 0, we can simply delete it from the X matrix. The new X matrix
# we'll call Xn, and it will be full rank.  We can use it to get betahat's
# and varbeta
Xn <- Xm[,-5]
Xn
betahat <- solve(t(Xn)%*%Xn)%*%t(Xn)%*%ym
betahat
sse <- t(ym-Xn%*%betahat)%*%(ym-Xn%*%betahat)
mse <- as.numeric(sse)/(15-5)
mse
varbeta <- mse*solve(t(Xn)%*%Xn)
varbeta
# While we cannot compute betahat from muhat since there are only 5 muhat's
# it turns out we can compute muhat from the betahat.  Since mu=ai%*%beta,
# we simply remove the column of ai associated with the degree of freedom
# we set equal to zero, and then multiply the result by betahat
newmuhat <- ai[,-5]%*%betahat
newmuhat
# and we have estimates for all 6 cell means
muhat
# we can also get a standard error for the muhat of the 1,1 cell
# using standard methods
vestmu <- ai[,-5]%*%varbeta%*%t(ai[,-5])
cbind(newmuhat,sqrt(diag(vestmu)))
# As we would expect the standard error for the first cell mean is quite a bit
# larger than that of the other cells since it was estimated using no data.
# Also, if we look at the linear combination represented by the first
# df for interaction, we will see that it is exactly 0
a[5,]%*%newmuhat
# We can also complete the F tests just as before, using our new muhat's
# First price
t(cp%*%newmuhat)%*%solve(cp%*%vestmu%*%t(cp))%*%cp%*%newmuhat/2
# advertising
t(ca%*%newmuhat)%*%solve(ca%*%vestmu%*%ca)%*%ca%*%newmuhat/1
# interaction
cin <- a[6,]
t(cin%*%newmuhat)%*%solve(cin%*%vestmu%*%cin)%*%cin%*%newmuhat/1
# only 1 df since we have constrained one to be exactly 0
# we are only testing the second interaction df
# Or you can use the x matrix that you created directly
t(betahat[2:3])%*%solve(varbeta[2:3,2:3])%*%betahat[2:3]/2
# test of adv
betahat[4]%*%solve(varbeta[4,4])%*%betahat[4]/1
# test of interaction
t(betahat[5])%*%solve(varbeta[5,5])%*%betahat[5]/1
