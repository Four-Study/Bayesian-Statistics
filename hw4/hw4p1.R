## load in necessary libraries
library(mvtnorm)

## Assume the data file is under the same directory as .R file
rat <- read.table("RatData.txt", row.names = 1)
colnames(rat) <- seq(1, 5)

########### Gibbs Sampler #############
# Gibbs.sampler <- function(data = NULL,
#                           a = 1,
#                           b = 1,
#                           eta = c(100, 40),
#                           Psi = matrix(c(1.1, -0.3, -0.3, 0.1), ncol = 2),
#                           rho = 30,
#                           R = matrix(c(1/6, 0.5, 0.5, 11/6), ncol = 2),
#                           verb = 0,
#                           burn.in = 1000) {
#   I <- nrow(data)
#   ## Initialize hyper-parameters
#   a <- 1
#   b <- 1
#   eta <- c(100, 40)
#   Psi <- matrix(c(1.1, -0.3, -0.3, 0.1), ncol = 2)
#   rho <- 30
#   R <- matrix(c(1/6, 0.5, 0.5, 11/6), ncol = 2)
#   
#   ## Initialize iteration 0
#   theta.init <- eta # modify for extreme case
#   phi.init <- a / (a + b) # modify for extreme case
#   Phi.init <- solve(R)
#   
#   ## MCMC process
#   M <- 1e6
#   theta0s <- rep(NA, M)
#   phis <- rep(NA, M)
#   Phis <- list
#   for (t in seq(1:M)) {
#     for (i in seq(1:I)) {
#       
#     }
#   }
# }

data <- rat
verb <- 10000
burn.in <- 1000

## Set hyper-parameters
a <- 0
b <- 0
eta <- c(100, 40)
Psi.inv <- matrix(0, ncol = 2, nrow = 2)
rho <- 2
R <- matrix(c(100, 0, 0, 0.1), ncol = 2)

## Initialize iteration 0
theta0.init <- eta # modify for extreme case
phi.init <- 0.1 # modify for extreme case
Phi.init <- matrix(c(1, 0, 0, 1), ncol = 2)

## Pre-process
I <- nrow(data)
J <- ncol(data)
X <- matrix(c(rep(1, J), seq(1:J)), ncol = 2)
sumX.1 <- t(X) %*% X

## MCMC process
M <- 1e5
theta0.s <- matrix(rep(NA, 2 * M), nrow = 2)
phi.s <- rep(NA, M)
Phi.s <- vector(mode = "list", length = M)
theta0.s[, 1] <- theta0.init
phi.s[1] <- phi.init
Phi.s[[1]] <- Phi.init

tik <- proc.time()
for (t in 2:M) {
  
  ## update theta_i's
  theta.i.s <- matrix(rep(NA, 2 * I), nrow = 2)
  for (i in 1:I) {
    sumX.2 <- t(as.matrix(data[i, ]) %*% X)
    sig <- solve(Phi.s[[t-1]] + phi.s[t-1] * sumX.1)
    mu <- solve(Phi.s[[t-1]] + phi.s[t-1] * sumX.1, 
                Phi.s[[t-1]] %*% theta0.s[, t-1] + phi.s[t-1] * sumX.2)
    theta.i.s[, i] <- rmvnorm(1, mean = mu, sigma = sig)
  }
  
  ## update theta_0
  sig <- solve(I * Phi.s[[t-1]] + Psi.inv)
  mu <- solve(I * Phi.s[[t-1]] + Psi.inv, 
              Phi.s[[t-1]] %*% rowSums(theta.i.s) + Psi.inv %*% eta)
  theta0.s[, t] <- rmvnorm(1, mean = mu, sigma = sig)
  
  ## update phi
  sumX.3 <- 0
  for (k in 1:I) 
    for (j in 1:J) 
      sumX.3 <- sumX.3 + (data[k, j] - t(X[j, ]) %*% theta.i.s[, i])^2
  sumX.3 <- sumX.3 / 2
  phi.s[t] <- rgamma(1, I*J/2 + a, sumX.3 + b)
  
  ## update Phi
  sig <- solve((theta.i.s - theta0.s[, t]) %*% 
                 t(theta.i.s - theta0.s[, t]) 
               + rho * R)
  Phi.s[[t]] <- matrix(rWishart(1, rho + I, sig), ncol = 2)
  
  ## monitor updating process
  if ((verb != 0) && (t %% verb ==0)) print(t)
}

tok <- proc.time()
cat(paste0("The traning process used ", (tok - tik)[3]), "\n")

## trace plots
plot(theta0.s[1, -(1:burn.in)], type = 'l')
plot(theta0.s[2, -(1:burn.in)], type = 'l')
plot(phi.s[-(1:burn.in)], type = 'l')

eigen.values <- matrix(rep(NA, 2 * M), nrow = 2)
traces <- rep(NA, M)
for (t in 1:M) {
  eigen.values[, t] <- eigen(Phi.s[[t]])$values
  traces[t] <- sum(diag(Phi.s[[t]]))
} 
plot(eigen.values[1, -(1:burn.in)], type = 'l')
plot(eigen.values[2, -(1:burn.in)], type = 'l')
plot(traces[-(1:burn.in)], type = 'l')
