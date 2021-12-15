######### Problem 2 #########

# Import the data and modify column names
dat <- read.csv("Model2_5444.txt", skip = 2, header = FALSE)
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

for (name in rownames(la.coef)[as.numeric(la.coef != 0) == 1]) {
  cat(name, "+ ")
}
la.coef[as.numeric(la.coef != 0) == 1]
which(as.numeric(la.coef != 0) == 1)
lm.fit1 <- lm(y ~ X3 + X13 + X17 + X27 + X45 + X2:X36 + 
                X7:X36 + X9:X14 + X9:X22 + X12:X44 + 
                X14:X19 + X15:X17 + X16:X35 + X18:X24 + 
                X20:X26 + X22:X46 + X23:X26 + X27:X44 + 
                X29:X32, data = dat.aug)
summary(lm.fit1)
coef.la <- lm.fit1$coefficients
library(tidyverse)
eqn.la <- paste("y =", paste(round(coef.la[1],3), 
                              paste(round(coef.la[-1],3), 
                                    names(coef.la[-1]), 
                                    sep=" * ", collapse=" + "), 
                              sep=" + ")) %>%
  gsub(pattern = ' \\* ', replacement = '*') %>%
  gsub(pattern = '\\+ -', replacement = '- ') %>%
  gsub(pattern = 'X([0-9]{1,2})', replacement = paste0("X_{", "\\1", "}")) %>%
  gsub(pattern = "[*:]", replacement = "")
print(eqn.la)

# AIC
null <- lm(y ~ 1, data = dat.aug)
lm.AIC <- step(null,scope=list(upper=terms(y ~ ., data=dat.aug)),
               trace = 0)
coef.AIC <- lm.AIC$coefficients
eqn.AIC <- paste("y =", paste(round(coef.AIC[1],3), 
                              paste(round(coef.AIC[-1],3), 
                                    names(coef.AIC[-1]), 
                                    sep=" * ", collapse=" + "), 
                              sep=" + ")) %>%
  gsub(pattern = ' \\* ', replacement = '*') %>%
  gsub(pattern = '\\+ -', replacement = '- ') %>%
  gsub(pattern = 'X([0-9]{1,2})', replacement = paste0("X_{", "\\1", "}")) %>%
  gsub(pattern = "\\*", replacement = "") %>%
  gsub(pattern = "\\.X", replacement = "X")
print(eqn.AIC)

# BIC
n <- nrow(dat.aug)
lm.BIC <- step(null,scope=list(upper=terms(y ~ ., data=dat.aug)),
               trace = 0, k = log(n))
coef.BIC <- lm.BIC$coefficients
eqn.BIC <- paste("y =", paste(round(coef.BIC[1],3), 
                              paste(round(coef.BIC[-1],3), 
                                    names(coef.BIC[-1]), 
                                    sep=" * ", collapse=" + "), 
                              sep=" + ")) %>%
  gsub(pattern = ' \\* ', replacement = '*') %>%
  gsub(pattern = '\\+ -', replacement = '- ') %>%
  gsub(pattern = 'X([0-9]{1,2})', replacement = paste0("X_{", "\\1", "}")) %>%
  gsub(pattern = "\\*", replacement = "") %>%
  gsub(pattern = "\\.X", replacement = "X")
print(eqn.BIC)

# SSVS

# initialization
p <- ncol(X)
psi <- 1
M <- 10000

Imat <- matrix(NA, nrow = M, ncol = p)
I <- c(TRUE, TRUE, rep(FALSE, p - 2))
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

plot(1:M, cumsum(t(Imat[, 50])) / 1:M, type = "l", ylim = c(0, 1))
