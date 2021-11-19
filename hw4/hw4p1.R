## load in necessary libraries
library(mvtnorm)

## Assume the data file is under the same directory as .R file
rat <- read.table("RatData.txt", row.names = 1)
colnames(rat) <- seq(1, 5)

## print iteration number every verb times
verb <- 1000

## Set up hyper-parameters
a <- 0
b <- 0
eta <- c(100, 50)
Psi.inv <- matrix(0, ncol = 2, nrow = 2)
rho <- 2
R <- matrix(c(100, 0, 0, 0.1), ncol = 2)

## Pre-process
I <- nrow(rat)
J <- ncol(rat)
X <- matrix(c(rep(1, J), seq(1:J)), ncol = 2)
sumX.1 <- t(X) %*% X

## MCMC process
M <- 1e4
theta0.s <- matrix(rep(NA, 2 * M), nrow = 2)
phi.s <- rep(NA, M)
Phi.s <- vector(mode = "list", length = M)

## Initialize iteration 1
theta0.s[, 1] <- eta
phi.s[1] <- 1
Phi.s[[1]] <- matrix(c(1, 0, 0, 1), ncol = 2)

## start iterations
for (t in 2:M) {
  
  ## update theta_i's
  theta.i.s <- matrix(rep(NA, 2 * I), nrow = 2)
  for (i in 1:I) {
    sumX.2 <- t(as.matrix(rat[i, ]) %*% X)
    sig <- solve(Phi.s[[t-1]] + phi.s[t-1] * sumX.1)
    mu <- sig %*% (Phi.s[[t-1]] %*% theta0.s[, t-1] + phi.s[t-1] * sumX.2)
    theta.i.s[, i] <- rmvnorm(1, mean = mu, sigma = sig)
  }
  
  ## update theta_0
  sig <- solve(I * Phi.s[[t-1]] + Psi.inv)
  mu <- sig %*% (Phi.s[[t-1]] %*% rowSums(theta.i.s) + Psi.inv %*% eta)
  theta0.s[, t] <- rmvnorm(1, mean = mu, sigma = sig)
  
  ## update phi
  sumX.3 <- 0
  for (i in 1:I) 
    for (j in 1:J) 
      sumX.3 <- sumX.3 + (rat[i, j] - t(X[j, ]) %*% theta.i.s[, i])^2
  sumX.3 <- sumX.3 / 2
  phi.s[t] <- rgamma(1, I*J/2 + a, sumX.3 + b)
  
  ## update Phi
  sig <- solve((theta.i.s - theta0.s[, t]) %*% 
                 t(theta.i.s - theta0.s[, t]) + rho * R)
  Phi.s[[t]] <- matrix(rWishart(1, rho + I, sig), ncol = 2)
  
  ## monitor updating process
  if ((verb != 0) && (t %% verb ==0)) print(t)
}

## trace plots
burn.in <- 100
pdf("imgs/theta0.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(4.1, 5.1, 1.1, 1.1))
plot(theta0.s[1, -(1:burn.in)], type = 'l', 
     xlab = "Iteration", ylab = "alpha 0", main = "", 
     cex.lab = 2, cex.axis = 1.5)
plot(density(theta0.s[1, -(1:burn.in)]), 
     ylab = "density", xlab = "alpha 0", main = "",
     cex.lab = 2, cex.axis = 1.5)
plot(theta0.s[2, -(1:burn.in)], type = 'l', 
     xlab = "Iteration", ylab = "beta 0", main = "",
     cex.lab = 2, cex.axis = 1.5)
plot(density(theta0.s[2, -(1:burn.in)]), 
     ylab = "density", xlab = "beta 0", main = "",
     cex.lab = 2, cex.axis = 1.5)
dev.off()

pdf("imgs/phi.pdf", width = 10, height = 5)
par(mar = c(4.1, 5.1, 1.1, 1.1))
plot(phi.s[-(1:burn.in)], type = 'l', 
     xlab = "Iteration", ylab = "phi", main = "", 
     cex.lab = 2, cex.axis = 1.5)
dev.off()

eigen.values <- matrix(rep(NA, 2 * M), nrow = 2)
traces <- Phi.s1 <- Phi.s2 <- rep(NA, M)
for (t in 1:M) {
  eigen.values[, t] <- eigen(Phi.s[[t]])$values
  traces[t] <- sum(diag(Phi.s[[t]]))
  Phi.s1[t] <- Phi.s[[t]][1, 1]
  Phi.s2[t] <- Phi.s[[t]][2, 2]
} 

pdf("imgs/Phi.pdf", width = 10, height = 10)
par(mfrow = c(2, 1), mar = c(4.1, 5.1, 1.1, 1.1))
plot(eigen.values[1, -(1:burn.in)], type = 'l', 
     xlab = "Iteration", ylab = "1st Eigenvalue", main = "", 
     cex.lab = 2, cex.axis = 1.5)
plot(eigen.values[2, -(1:burn.in)], type = 'l', 
     xlab = "Iteration", ylab = "2nd Eigenvalue", main = "", 
     cex.lab = 2, cex.axis = 1.5)
dev.off()

save.image("p1.RData")
