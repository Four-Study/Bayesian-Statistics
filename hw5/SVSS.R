# Import the data and modify column names
dat <- read.csv("Model1_5444.txt", skip = 2, header = FALSE)
colnames(dat) <- c("y", paste0(rep("X", 50), 1:50))

# augment the data set with interactions
X <- model.matrix(y ~ .^2, dat)
y <- dat[, 1]
dat.aug <- data.frame(cbind(y = y, X[, -1]))

# LASSO
library(glmnet)
lambda <- exp(seq(-10, 1, 0.1))
la.fit1 <- glmnet(X, y, alpha = 1, lambda = lambda)
cv.out <- cv.glmnet(X, y, alpha = 1, lambda = lambda)
best.lam <- cv.out$lambda.min
la.coef <- predict(la.fit1,type="coefficients",s=best.lam)
la.coef[as.numeric(la.coef != 0) == 1]
init <- which(as.numeric(la.coef != 0) == 1)

# SSVS

# initialization
p <- ncol(X)
psi <- 1
M <- 10000

Imat <- matrix(NA, nrow = M, ncol = p)
I <- rep(FALSE, p)
I[init] <- TRUE
for (i in 1:M) {
  for (j in 1:p) {
    # calculate the numerator and denominator of Bayes factor
    I.tmp <- I
    I.tmp[j] <- FALSE
    X.r <- X[, I.tmp]
    inv.r <- solve(t(X.r) %*% X.r + 1 / psi^2 * diag(max(ncol(X.r), 1)))
    log.reduced <- log(psi) + y %*% X.r %*% inv.r %*% t(X.r) %*% y +  
      log(sqrt(det(inv.r)))
    I.tmp[j] <- TRUE
    X.f <- X[, I.tmp]
    inv.f <- solve(t(X.f) %*% X.f + 1 / psi^2 * diag(ncol(X.f)))
    log.full <- y %*% X.f %*% inv.f %*% t(X.f) %*% y + log(sqrt(det(inv.f)))
    
    # calculate Bayes factor
    BF <- exp(log.reduced - log.full)
    prob <- 1 / (1 + 1 / BF)
    I[j] <- rbinom(1, 1, prob) == 0
  }
  Imat[i, ] <- I
  if ((i %% 100) == 0) cat("loop", i, "\n")
}
# plot(1:M, cumsum(t(Imat[, 11])) / 1:M, type = "l", ylim = c(0, 1))
save.image("1.RData")
