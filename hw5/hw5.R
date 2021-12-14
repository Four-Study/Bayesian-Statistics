######### Problem 2 #########

# Import the data and modify column names
dat <- read.csv("Model1_5444.txt", skip = 2, header = FALSE)
colnames(dat) <- c("Y", paste0(rep("X", 50), 1:50))

# augment the data set with interactions
all.x <- model.matrix(Y ~ .^2, dat)
dat.aug <- data.frame(cbind(Y = dat[, 1], all.x[, -1]))

# LASSO
library(glmnet)
lambda <- exp(seq(-10, 1, 0.1))
la.fit1 <- glmnet(dat.aug[, -1], dat.aug[, 1], alpha = 1, lambda = lambda)


plot(la.fit1$lambda, la.fit1$dev.ratio, type = 'l')
plot(la.fit1)
s <- 0.3
rownames(coef(la.fit1, s = s))[as.numeric(coef(la.fit1, s = s) != 0) == 1]
lm.fit1 <- lm(Y ~ X3 + X17 + X27 + X45 + 
                X3.X15 + X31.X45, data = dat.aug)
summary(lm.fit1)

# AIC
library(MASS)
null <- lm(Y ~ 1, data = dat[-(49:50)])
full <- lm(Y ~ ., data = dat[-(49:50)]) # dat.aug
stepAIC(object = full, scope = (~ .), trace = FALSE)
AIC(null)
AIC(full)

# BIC
n <- nrow(dat.aug)
stepAIC(object = full, scope = ( ~ .), k = log(n),
        trace = 1)

# SSVS

# initialization
X <- as.matrix(cbind(1, dat[, -(1:2)]))
y <- dat[, 1]
p <- ncol(X)
psi <- 1
M <- 10000

Imat <- matrix(NA, nrow = M, ncol = p)
I <- rep(TRUE, p)
for (i in 1:M) {
  for (j in 1:p) {
    # calculate the numerator and denominator of Bayes factor
    I.tmp <- I
    I.tmp[j] <- FALSE
    X.r <- X[, I.tmp]
    inv.r <- solve(t(X.r) %*% X.r + 1 / psi^2 * diag(ncol(X.r)))
    log.reduced <- log(psi) + y %*% X.r %*% inv.r %*% t(X.r) %*% y +  
      log(sqrt(det(inv.r)))
    I.tmp[j] <- TRUE
    X.f <- X[, I.tmp]
    inv.f <- solve(t(X.f) %*% X.f + 1 / psi^2 * diag(ncol(X.f)))
    log.full <- y %*% X.f %*% inv.f %*% t(X.f) %*% y + log(sqrt(det(inv.f)))
    
    # calculate Bayes factor
    if (I[j]) BF <- exp(log.reduced - log.full)
    else BF <- exp(log.full - log.reduced)
    prob <- 1 / (1 + 1 / BF)
    I[j] <- rbinom(1, 1, prob) == 0
  }
  Imat[i, ] <- I
  if ((i %% 100) == 0) cat("loop", i, "\n")
}

plot(1:M, cumsum(t(Imat[, 50])) / 1:M, type = "l", ylim = c(0, 1))
