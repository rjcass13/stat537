library('nlme')
library('lme4')

weight <- read.table('weight.dat', header = TRUE)

# Tmt: 
# CONT: Control: 20 subjects
# RI: Rep increase: 16 subjects
# Weight increase: 21 subjects

# Time: 7 time intervals, each subject measured at each interval

# Make each subject unique
weight$sub <- paste(substr(weight$tmt, 1, 1), sprintf("%02d", weight$subject), sep='')
# Make time start at 0. Lots of reasons, but for one if they don't start at the same spot, 
# there was probably an issue with randomization
weight$ztime <- weight$time - 1
weightsort <- weight[order(weight$sub), ]
weightsort[1:14, ]

# Some nice quick views of the data. Look at intercepts, slopes, etc. Note C14 is heavy negative slope
WeightG <- groupedData(response ~ time | sub, data = weight)
plot(WeightG)
plot(WeightG, outer = ~tmt)

# This model gives 3 intercepts and 3 slopes (we get rid of the normal intercept to get values on ALL groups)
# This structures the design matrix for us to actually see the slopes of all groups
fitlm <- lm(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt), data = weight)
# 4th line indicates control doesn't really change over time
summary(fitlm)
# However, this is the wrong analysis. It's not accounting for subject at all, etc.
logLik(fitlm)

# Fit mixed model. By default fits CompoundSymmetric
# Does CompSymm make sense for this data? It implies that with a subject, the covariance between 
# times 1&2 is the same as 1&7, etc.
fitcs <- lmer(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt) + (1|sub), data = weight)
summary(fitcs)
logLik(fitcs)
# We're fitting 3 slopes, but 57 lines (the lines just have slope according to the treatment, but intercept 
# according to the subject)
# X is a matrix that indicates: first 3 columns which treatment per time step,
#  ext 3 columns the time steps per treatment)
# beta is the intercept
# Z is matrix that indicates which treatment the subjects were in
# u is the deviation from the intercept
# This structure allows each line (57 subjects) to start at a different spot

# Same exact structure as the model above, just different function
fitcsg <- gls(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt), 
  correlation = corCompSymm(form = ~1|sub), data = weight)
logLik(fitcsg)
# However, fitcs gives us the u-hats, fitcsg does not. If we care about ind. subjects, we might want u-hats
# You can see the u-hats using:
ranef(fitcs)

# The random effect is now putting both the intercept and time in the Z-matrix
# This gives us a deviation from the intercept AND a deviation from the slope
# Now each line can start at a different place and tilt differently
fitrc <- lmer(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt) + (ztime|sub), data = weight)
summary(fitrc)
# u ~ N(0,G). G was a diagonal with sigma^2_s
# Now in the summary we see there is a Corr of -0.06 on ztime
# So now u has off-diagonals??
logLik(fitrc)

# The || assumes the slopes and intercepts are independent (not correlation)
# Generally we'd assume slopes and intercepts have negative correlation (slope up, intercept down)
fitrc1 <- lmer(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt) + (ztime||sub), data = weight)
summary(fitrc1)
logLik(fitrc1)

# ANOVA requires ML, not REML. In REML we integrate out the betas, but the ANOVA needs it
anova(fitrc, fitrc1)
logLik(fitrc) #-659.6185 (df = 10)
logLik(fitrc1) # -659.6982 (df = 9), diff ~.08, *2 ~ .16
1-pchisq(.16, 1) # Not significantly differently


bhat <- t(t(fixef(fitrc1)))
uhat <- ranef(fitrc1)$sub # 2 columns, intercepts and slopes

# Using y = XB + Zu + e, e ~ N(0, R), u ~ N(0, G), e and u are independent
# We will use the C-hat shown in the SAS documentation
uhat <- matrix(as.matrix(uhat), 114, 1) # combines it all into one column
X <- model.matrix(~ -1 + as.factor(tmt) + ztime:as.factor(tmt), data = weightsort)
dim(X) # 399 x 6
Z <- model.matrix(response ~ -1 + as.factor(sub) + ztime:sub, data = weightsort)
dim(Z) # 399 x 114
Z[2, ] # Start with a 1 for intercept of C01, and a 1 for the slope of C01
summary(fitrc1)
R <- .6337*diag(1,399) # Residual variance x Identity
g <- rbind(c(8.8853, 0), c(0, 0.1234)) # Intercept then Slope
g # This is our little 2x2 block
G <- g %x% diag(1, 57)
G[1:7, 1:7] # First 57 rows/columns is a diagonal of 8.8853
G[56:60, 56:60] # Last 57 rows/columns is a diagonal of 0.1234. This shows the transition

sigma <- Z %*% G %*% t(Z) + R
sigma[1:7, 1:7]

sigma_i <- solve(sigma)
v11 <- solve(t(X) %*% sigma_i %*% X)
v22 <- solve(t(Z) %*% solve(R) %*% Z + solve(G)) + G %*% t(Z) %*% sigma_i %*% X %*% v11 %*% t(X) %*% sigma_i %*% Z %*% G
v12 <- -v11 %*% t(X) %*% sigma_i %*% Z %*% G
v21 <- t(v12)
chat <- rbind(cbind(v11, v12), cbind(v21,v22))

z56 <- t(rep(0,56))

# This is getting the value for first subject on 6th (last) time period
# 3 intercept columns, 3 slope columns
# Should be 79.97 - .048*6 + 4.95 + .344*6, pulled from the summary
csub1 <- cbind(1,0,0,6,0,0,1,z56,6,z56)
bhatuhat <- rbind(bhat, uhat)
csub1 %*% bhatuhat

# Intervals for subject 1 in control at ztime 6
fit_test <- lmerTest::lmer(response ~ as.factor(tmt) + ztime:as.factor(tmt) + (ztime||sub), data = weightsort)
anova_table <- anova(fit_test)
den_df <- anova_table$DenDF[1]  # Denominator degrees of freedom

# Confidence Interval
csub1 %*% bhatuhat + qt(.975, den_df)*sqrt(csub1 %*% chat %*% t(csub1))
csub1 %*% bhatuhat - qt(.975, den_df)*sqrt(csub1 %*% chat %*% t(csub1))

# Prediction interval (which is a little silly here)
# Uses the residual variance (e)
csub1 %*% bhatuhat + qt(.975, den_df)*sqrt(csub1 %*% chat %*% t(csub1) + .6337)
csub1 %*% bhatuhat - qt(.975, den_df)*sqrt(csub1 %*% chat %*% t(csub1) + .6337)

# Now looking at the overal group estimate for treatment 3 at time 6
vhatbhat <- chat[1:6,1:6]
cw6 <- rbind(c(0,0,1,0,0,6))
ucl <- cw6 %*% bhat + qt(.975, den_df)*sqrt(cw6 %*% vhatbhat %*% t(cw6))
lcl <- cw6 %*% bhat - qt(.975, den_df)*sqrt(cw6 %*% vhatbhat %*% t(cw6))
c(lcl, cw6%*%bhat, ucl)

# Now looking at all times for treatment 3
cw <- cbind(0,0,1,0,0,0:6)
plot(1:7, cw%*%bhat, xlim = c(0,8), ylim = c(75,90), type = "l")
ucl <- cw %*% bhat + qt(.975, den_df)*sqrt(diag(cw %*% vhatbhat %*% t(cw)))
lcl <- cw %*% bhat - qt(.975, den_df)*sqrt(diag(cw %*% vhatbhat %*% t(cw)))
cbind(lcl, cw%*%bhat, ucl)
lines(1:7, ucl, col = "green")
lines(1:7, lcl, col = "green")

# Now add on the plots for treatment 2
cr <- cbind(0,1,0,0,0:6,0)
lines(1:7, cr%*%bhat, lty = 3, col = "blue")
ucl_r <- cr %*% bhat + qt(.975, den_df)*sqrt(diag(cr %*% vhatbhat %*% t(cr)))
lcl_r <- cr %*% bhat - qt(.975, den_df)*sqrt(diag(cr %*% vhatbhat %*% t(cr)))
cbind(lcl_r, cr%*%bhat, ucl_r)
lines(1:7, ucl_r, col = "orange", lty = 3)
lines(1:7, lcl_r, col = "orange", lty = 3)

# Test the difference in the intercepts
crw1 <- cbind(0,-1,1,0,0,0)
tstat <- (crw1 %*% bhat)/sqrt(crw1 %*% vhatbhat %*% t(crw1))
2 * (1-pt(tstat, den_df))

# Test the difference in the slopes
csrw <- cbind(0,0,0,0,-1,1)
tstat <- (csrw %*% bhat)/sqrt(csrw %*% vhatbhat %*% t(csrw))
2 * (1-pt(tstat, den_df))

# Now the prediction interval lines
upl_r <- cr %*% bhat + qt(.975, den_df)*sqrt(diag(cr %*% vhatbhat %*% t(cr) + .6337))
lpl_r <- cr %*% bhat - qt(.975, den_df)*sqrt(diag(cr %*% vhatbhat %*% t(cr) + .6337))
lines(1:7, upl_r, col = "green")
lines(1:7, lpl_r, col = "green")

upl_w <- cw %*% bhat + qt(.975, den_df)*sqrt(diag(cw %*% vhatbhat %*% t(cw) + .6337))
lpl_w <- cw %*% bhat - qt(.975, den_df)*sqrt(diag(cw %*% vhatbhat %*% t(cw) + .6337))
lines(1:7, upl_w, col = "orange")
lines(1:7, lpl_w, col = "orange")

# Now use the variance from the data
z <- cbind(1,0:6)
r <- .6337*diag(1,7)
vardata <- z %*% g %*% t(z) + r
upl_rr <- cr %*% bhat + qt(.975, den_df) * sqrt(diag(cr %*% vhatbhat %*% t(cr) + vardata))
lpl_rr <- cr %*% bhat + qt(.975, den_df) * sqrt(diag(cr %*% vhatbhat %*% t(cr) + vardata))
lines(1:7, upl_rr, col = "orange", lty = 3)
lines(1:7, lpl_rr, col = "orange", lty = 3)

# Other potential models
fitar1 <- gls(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt), correlation = corAR1(form = ~ ztime | sub), data = weightsort)
summary(fitar1)
3.282^2 # Variance
logLik(fitar1)
getVarCov(fitar1)
vhatbhat <- vcov(fitar1)
vhatbhat
bb <- matrix(0,7,7)

for (i in 1:7) { 
  for (j in 1:7) {
    bb[i,j] <- ifelse(i==j, 10.77, 10.77*.95^(abs(i-j)))
  }
}
fitar1w <- gls(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt), correlation = corAR1(form = ~ ztime | sub), weight = varIdent(form = ~1|tmt), data = weight)
# This makes each block covariance (for each treatment) has it's own diagonal (scaled sigma2_e)
summary(fitar1w) 
logLik(fitar1w)
anova(fitar1, fitar1w)
# Better fit, but we would assume WI is slightly bigger than RI (higher weights should have higher variance),
# which is not the case here. Doesn't pass the "sniff test"
fitar1wzt <- gls(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt), correlation = corAR1(form = ~ ztime | sub), weight = varIdent(form = ~1|ztime), data = weight)
summary(fitar1wzt)
anova(fitar1, fitar1wzt)
logLik(fitar1wzt)
getVarCov(fitar1wzt)
# Unstructured
fitun <- gls(response ~ -1 + as.factor(tmt) + ztime:as.factor(tmt), correlation = corSymm(form = ~1|sub), weight = varIdent(form = ~1|ztime), data = weight)
logLik(fitun)
getVarCov(fitun)
anova(fitar1, fitun)
summary(fitun)

AIC(fitrc1)
AIC(fitar1)
AIC(fitun)
AIC(fitar1w)
AIC(fitar1wzt)
BIC(fitrc1)
BIC(fitar1)

keep <- sample(1:399, 300)
weightun <- weight[keep,]
dim(weight)
fitrc1u <- lmer(response ~ -1 + as.factor(tmt)+ztime:as.factor(tmt) + (ztime||sub),data=weightun)
fitar1u <- gls(response ~-1+as.factor(tmt)+ztime:as.factor(tmt),correlation = corAR1(form=~1|sub),data=weightun)
fitar1wu <- gls(response ~-1+as.factor(tmt)+ztime:as.factor(tmt),correlation = corAR1(form=~1|sub),weights=varIdent(form=~1|tmt),data=weightun)
AIC(fitrc1u)
AIC(fitar1u)
AIC(fitar1wu)
BIC(fitrc1u)
BIC(fitar1u)
BIC(fitar1wu)