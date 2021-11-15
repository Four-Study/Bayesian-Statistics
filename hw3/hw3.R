## Problem 1
require(mvtnorm)
n <- 1000
mu <- 0
Sig <- matrix(c(1, 0.8, 0.8, 1), ncol = 2)
nu <- 1
set.seed(12345)
xy <- rmvt(n, sigma=Sig, df=nu)
par(mar=c(5,5,4,2))
plot(xy, xlab = "x", ylab = "y", cex.axis = 1.5, cex.lab = 2)

lm_fit <- lm(xy[, 2] ~ xy[, 1])
par(mfrow = c(1, 2))
plot(lm_fit$fitted.values, lm_fit$residuals, 
     xlab = "Fitted Values",
     ylab = "Residuals", 
     cex.axis = 1.5, cex.lab = 2)
qqnorm(lm_fit$residuals, cex.axis = 1.5, cex.lab = 2, main = "")
qqline(lm_fit$residuals, col = "red")
par(mfrow = c(1, 1))
