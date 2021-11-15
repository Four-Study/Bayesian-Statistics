################# part a & b ###############
# generate the probability sequence
ps <- seq(0.01, 0.99, 0.01)
# define the function for sampling
sim <- function(p, N = 30) {
  xs <- rbinom(1000, N, p)
  return(xs)
}
## Note: change N to 50, 100 and 1000 for part b
N <- 30
xss <- sapply(ps, function(x){sim(x, N = N)})
phat <- xss / N
# lower and upper bounds
conf_lower <- phat - 1.96 * sqrt(phat * (1 - phat) / N)
conf_upper <- phat + 1.96 * sqrt(phat * (1 - phat) / N)
# if the interval covers the truth
conf_cover <- matrix(NA, 1000, length(ps))
for (i in 1:length(ps)) {
  conf_cover[, i] <- (conf_upper[, i] >= ps[i]) & (conf_lower[, i] <= ps[i])
}
conf_freq <- apply(conf_cover, 2, sum)

# The prior distribution is Beta(1/2, 1/2)
cred_lower <- qbeta(0.025, xss + 1/2, N - xss + 1/2)
cred_upper <- qbeta(0.975, xss + 1/2, N - xss + 1/2)
# if the interval covers the true value
cred_cover <- matrix(NA, 1000, length(ps))
for (i in 1:length(ps)) {
  cred_cover[, i] <- (cred_upper[, i] >= ps[i]) & (cred_lower[, i] <= ps[i])
}
cred_freq <- apply(cred_cover, 2, sum)

par(mar = c(5, 5, 5, 5))
plot(ps, conf_freq, type = 'l', lwd = 2, 
     ylim = c(200, 1000), 
     xlab = substitute(paste("Value of ", italic('p'))),
     ylab = "Coverage Frequency",
     cex.axis = 1.5,
     cex.lab = 2)
lines(ps, cred_freq, lwd = 2, lty = 2)
legend("bottom", 
       c("Confidence Interval Coverage Frequency", 
         "Credible Interval Coverage Frequency"), 
       lty = c(1, 2), cex = 1.6)

############### part c & d ##################
set.seed(12345)
p <- 0.001
Ns <- seq(5, 4000, 5)
xss2 <- sapply(Ns, function(x){sim(p = p, N = x)})

phat2 <- matrix(NA, 1000, length(Ns))
for (i in 1:length(Ns)) {
  phat2[, i] <- xss2[, i] / Ns[i]
}
# lower and upper bounds
conf_lower2 <- phat2 - 1.96 * sqrt(phat2 * (1 - phat2) / N)
conf_upper2 <- phat2 + 1.96 * sqrt(phat2 * (1 - phat2) / N)
# if the interval covers the truth
conf_cover2 <- (conf_upper2 >= p) & (conf_lower2 <= p)
conf_rate2 <- apply(conf_cover2, 2, sum) / 1000

# The prior distribution is Beta(1/2, 1/2)
cred_lower2 <- matrix(NA, 1000, length(Ns))
cred_upper2 <- matrix(NA, 1000, length(Ns))
for (i in 1:length(Ns)) {
  cred_lower2[, i] <- qbeta(0.025, xss2[, i] + 1/2, Ns[i] - xss2[, i] + 1/2)
  cred_upper2[, i] <- qbeta(0.975, xss2[, i] + 1/2, Ns[i] - xss2[, i] + 1/2)
}

# if the interval covers the true value
cred_cover2 <- (cred_upper2 >= p) & (cred_lower2 <= p)
cred_rate2 <- apply(cred_cover2, 2, sum) / 1000
# plot
plot(Ns, conf_rate2, type = 'l', lwd = 2, 
     xlab = substitute(paste("Value of ", italic('N'))),
     ylab = "Coverage Rate",
     cex.axis = 1.5,
     cex.lab = 2)
lines(Ns, cred_rate2, lwd = 2, lty = 7, col = "gray")
abline(a = 0.95, b = 0, col = "red")
legend("bottom", 
       c("Confidence Interval Coverage Rate", 
         "Credible Interval Coverage Rate",
         "0.95 Nominal Rate"), 
       col = c("black", "gray", "red"),
       lty = c(1, 2), cex = 1.6)
# Calculation
idx1 <- min(which(abs(conf_rate2 - 0.95) < 0.01))
Ns[idx1]
idx2 <- min(which(abs(cred_rate2 - 0.95) < 0.01))
Ns[idx2]
