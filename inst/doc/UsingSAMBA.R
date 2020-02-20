## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, fig.cap="Figure 1: Model Structure", out.width = '80%'------
knitr::include_graphics("images/ModelDiagram.png")

## ---- echo = TRUE, eval = TRUE,  fig.width = 7, fig.height = 4----------------

library(SAMBA)
library(MASS)
expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x / (1 - x))

nobs <- 5000

### Generate Predictors and Follow-up Information
set.seed(1234)
cov <- mvrnorm(n = nobs, mu = rep(0, 3), Sigma = rbind(c(1,   0, 0.4),
                                                       c(0,   1,   0),
                                                       c(0.4, 0,   1)))

data <- data.frame(Z = cov[, 1], X = cov[, 2], W = cov[, 3])

# Generate random uniforms
set.seed(5678)
U1 <- runif(nobs)
set.seed(4321)
U2 <- runif(nobs)
set.seed(8765)
U3 <- runif(nobs)

# Generate Disease Status
DISEASE <- expit(-2 + 0.5 * data$Z)
data$D   <- ifelse(DISEASE > U1, 1, 0)

# Relate W and D
data$W <- data$W + 1 * data$D

# Generate Misclassification
SENS <- expit(-0.4 + 1 * data$X)
SENS[data$D == 0] = 0
data$Dstar <- ifelse(SENS > U2, 1, 0)

# Generate Sampling Status
SELECT <- expit(-0.6 + 1 * data$D + 0.5 * data$W)
S  <- ifelse(SELECT > U3, T, F)

# Observed Data
data.samp <- data[S,]

# True marginal sampling ratio
prob1 <- expit(-0.6 + 1 * 1 + 0.5 * data$W)
prob0 <- expit(-0.6 + 1 * 0 + 0.5 * data$W)
r.marg.true <- mean(prob1[data$D == 1]) / mean(prob0[data$D == 0])

# True inverse probability of sampling weights
prob.WD <- expit(-0.6 + 1 * data.samp$D + 0.5 * data.samp$W)
weights <- nrow(data.samp) * (1  / prob.WD) / (sum(1 / prob.WD))

# True associations with D in population
trueX <- glm(D ~ X, binomial(), data = data)
trueZ <- glm(D ~ Z, binomial(), data = data)

# Initial Parameter Values
fitBeta  <- glm(Dstar ~ X, binomial(), data = data.samp)
fitTheta <- glm(Dstar ~ Z, binomial(), data = data.samp)

## ---- results='hide', message=F-----------------------------------------------
# Using marginal sampling ratio r and P(D=1)
sens1 <- sensitivity(data.samp$Dstar, data.samp$X, mean(data$D),
                      r = r.marg.true)

# Using inverse probability of selection weights and P(D=1)
sens2 <- sensitivity(data.samp$Dstar, data.samp$X, prev = mean(data$D),
                     weights = weights)

# Using marginal sampling ratio r and P(D=1|X)
prev  <- predict(trueX, newdata = data.samp, type = 'response')
sens3 <- sensitivity(data.samp$Dstar, data.samp$X, prev, r = r.marg.true)

# Using inverse probability of selection weights and P(D=1|X)
prev  <- predict(trueX, newdata = data.samp, type = 'response')
sens4 <- sensitivity(data.samp$Dstar, data.samp$X, prev, weights = weights)

## ---- results='hide'----------------------------------------------------------
# Approximation of D*|Z
approx1 <- approxdist(data.samp$Dstar, data.samp$Z, sens1$c_marg,
                      weights = weights)

# Non-logistic link function method
nonlog1 <- nonlogistic(data.samp$Dstar, data.samp$Z, c_X = sens3$c_X,
                       weights = weights)

# Direct observed data likelihood maximization without fixed intercept
start <- c(coef(fitTheta), logit(sens1$c_marg), coef(fitBeta)[2])
fit1 <- obsloglik(data.samp$Dstar, data.samp$Z, data.samp$X, start = start,
                 weights = weights)
obsloglik1 <- list(param = fit1$param, variance = diag(fit1$variance))

# Direct observed data likelihood maximization with fixed intercept
fit2   <- obsloglik(data.samp$Dstar, data.samp$Z, data.samp$X, start = start,
                 beta0_fixed = logit(sens1$c_marg), weights = weights)
obsloglik2 <- list(param = fit2$param, variance = diag(fit2$variance))

# Expectation-maximization algorithm without fixed intercept
fit3 <- obsloglikEM(data.samp$Dstar, data.samp$Z, data.samp$X, start = start,
                 weights = weights)
obsloglik3 <- list(param = fit3$param, variance = diag(fit3$variance))

# Expectation-maximization algorithm with fixed intercept
fit4 <- obsloglikEM(data.samp$Dstar, data.samp$Z, data.samp$X, start = start,
                  beta0_fixed = logit(sens1$c_marg), weights = weights)
obsloglik4 <- list(param = fit4$param, variance = diag(fit4$variance))

## ---- echo = FALSE, eval = TRUE,  fig.width = 5, fig.height= 5----------------
plot(sort(sens3$c_X), xlab = 'Patients', ylab = 'Sensitivity',
     main = 'Figure 2: Sensitivity Estimates', type = 'l', col = 'red', lwd = 2)
lines(sort(expit(obsloglik1$param[3] + obsloglik1$param[4]*data.samp$X)), col = 'blue', lwd = 2)
lines(sort(expit(obsloglik2$param[3] + obsloglik2$param[4]*data.samp$X)), col = 'green', lwd = 2)
abline(h=sens1$c_marg, col = 'purple', lwd = 2)
lines(sort(expit(-0.4 + 1*data.samp$X)), col = 'black', lwd = 2)
legend(x='topleft', fill = c('purple', 'red','blue', 'green', 'black'),
       legend = c('Estimated marginal sensitivity',
                  'Using non-logistic link method',
                  'Using obs. data log-lik',
                  'Using obs. data log-lik (fixed intercept)',
                  'Truth'), cex = 0.7)

## ---- echo = FALSE, eval = TRUE,  fig.width = 5, fig.height= 5,  results='hide', message=F----
rvals = c(1,1.5,2,2.5,5,10)
COL = c('red', 'orange', 'yellow', 'green', 'blue', 'purple')
true_prevs = predict(trueX, newdata = data.samp, type = 'response')
plot(sort(expit(-0.4 + 1*data.samp$X)), xlab = 'Patients', ylab = 'Sensitivity',
main = 'Figure 3: Estimated sensitivity across \n marginal sampling ratios',
type = 'l', col = 'black', lwd = 2, ylim = c(0,1))

for (i in 1:length(rvals)) {
  TEMP <- sensitivity(X = data.samp$X, Dstar = data.samp$Dstar,  r = rvals[i], prev = true_prevs)
  lines(sort(TEMP$c_X), col = COL[i])
}
legend(x='topleft', legend = c(rvals, 'Truth'), title = 'Sampling Ratio',
       fill = c(COL, 'black'), cex = 0.8)

## ---- echo = FALSE, eval = TRUE,  fig.width = 5, fig.height= 5----------------
plot(fit1$beta0_fixed, fit1$loglik.seq, xlab = 'Beta_0', ylab = 'Log-likelihood',
     main = 'Figure 4: Profile Log-Likelihood Values \n for Direct Maximization', pch = 16)
abline(v=-0.4, col = 'blue', lwd = 2)
points(logit(sens1$c_marg), fit2$loglik.seq, pch = 17, col = 'red')
legend(x='topright', legend = c('No fixed beta_0', 'Fixed beta_0'), col = c('black', 'red'), pch = c(16,17))
text(x=-0.1, y=mean(fit1$loglik.seq),label = 'True beta_0', srt = 90)

plot(fit3$loglik.seq[-c(1)], xlab = 'EM algorithm iteration', ylab = 'Log-likelihood',
     main = 'Figure 5: Log-Likelihood Values \n Across EM Iterations', pch = 16)
points(fit4$loglik.seq[-c(1)], pch = 17, col = 'red')
legend(x='bottomright', legend = c('No fixed beta_0', 'Fixed beta_0'), pch = c(16,17),
       col = c('black', 'red'))

## ---- echo = FALSE, eval = TRUE,  fig.width = 7, fig.height= 4, message=FALSE----
library(ggplot2)
library(scales)

## ---- echo = FALSE, eval = TRUE,  fig.width = 7, fig.height= 4----------------
METHODS = c('True',  'Uncorrected', 'Approx D*|Z + IPW','Non-logistic Link + IPW','Obs. log-lik + IPW', 'Fixed intercept obs. log-lik + IPW')
PARAM = c( coef(trueZ)[2], coef(fitTheta)[2],approx1$param,  nonlog1$param[2], obsloglik1$param[2], obsloglik2$param[2] )
VARIANCE = c(diag(summary(trueZ)$cov.scaled)[2],diag(summary(fitTheta)$cov.scaled)[2],
             approx1$variance,nonlog1$variance[2],obsloglik1$variance[2], obsloglik2$variance[2])
pd = position_dodge(width=0.6)
a <- ggplot(data = data.frame(METHODS = METHODS, PARAM = PARAM, VARIANCE = VARIANCE),
       aes(xmin= METHODS, xmax = METHODS, ymin = PARAM - 1.96*sqrt(VARIANCE), ymax =  PARAM + 1.96*sqrt(VARIANCE),
           col = METHODS,x = METHODS, y = PARAM)) +
  geom_point(position = position_dodge(.7), size = 2) +
  geom_linerange(position = position_dodge(.7), size = 1.2) +
  xlab('') + ylab('logOR')+ggtitle('Figure 6: Estimated Log-Odds Ratio Across Methods')+
  scale_x_discrete(limits=METHODS)+
  geom_hline(yintercept = PARAM[1], linetype = 1, color = 'black')+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  theme(legend.position="top",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.text=element_text(size=8), legend.title = element_blank(),
        axis.text.x=element_text(angle=20,hjust=1,vjust=1), text = element_text(size=12))
print(a)

