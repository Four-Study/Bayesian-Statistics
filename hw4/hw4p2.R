library(mvtnorm)

N <- 1000
mu <- 0
Sig <- matrix(c(1, 0.8, 0.8, 1), ncol = 2)
nu <- 1
set.seed(12345)
xy <- rmvt(N, sigma=Sig, df=nu)
x <- xy[, 1]
X <- matrix(c(rep(1, N), x), ncol = 2)
y <- xy[, 2]

verb <- 10000
burn.in <- 1000

## set hyper-parameters
a <- 0.5
b <- 0.5

## MCMC process
M <- 1e5
beta.s <- matrix(rep(NA, 2 * M), nrow = 2)
phi.s <- rep(NA, M)
beta.s[, 1] <- c(1, 2)
phi.s[1] <- 1

tik <- proc.time()

for (t in 2:M) {
  
  ## update gamma_i's
  resids <- y - beta.s[1, t-1] - beta.s[2, t-1]*x
  gamma.i.s <- rgamma(N, a + 0.5, 1/2 * phi.s[t-1] * resids^2 + b)
  
  ## update beta
  Gamma <- diag(gamma.i.s)
  sig <- solve(phi.s[t-1] * t(X) %*% Gamma %*% X)
  mu <- solve(t(X) %*% Gamma %*% X, t(X) %*% Gamma %*% y)
  beta.s[, t] <- rmvnorm(1, mu, sig)
  
  ## update phi
  resids <- y - X %*% beta.s[, t]
  rate <- 1/2 * t(resids) %*% Gamma %*% resids
  phi.s[t] <- rgamma(1, N / 2, rate)
  
  ## monitor iterating process
  print(t)
  if ((verb != 0) && (t %% verb ==0)) print(t)
}

tok <- proc.time()
cat(paste0("The traning process used ", (tok - tik)[3]), "\n")

## Visualization
plot(beta.s[1, 1:t], type = 'l')
plot(beta.s[2, 1:t], type = 'l')
