## problem 2
mu <- 5
sig <- 1
N <- 1000
result <- matrix(NA, nrow = 5, ncol = 2)
X <- rep(NA, N)
set.seed(20210928)
for (k in 1:5) {
  for (i in 1:N) {
    gamma <- rgamma(1, 1/2, rate = 1/2)
    X[i] <- rnorm(1, mu, sig / sqrt(gamma))
  }
  result[k, 1] <- min(X)
  result[k, 2] <-max(X)
}
