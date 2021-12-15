# Import the data and modify column names
dat <- read.csv("Model1_5444.txt", skip = 2, header = FALSE)
colnames(dat) <- c("y", paste0(rep("X", 50), 1:50))

# SSVS

SVSS <- function(X, y, M = 10000) {
  # initialization
  p <- ncol(X)
  psi <- 1
  
  Imat <- matrix(NA, nrow = M, ncol = p)
  I <- rep(TRUE, p)
  # I[init] <- TRUE
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
    if ((i %% 1000) == 0) cat("loop", i, "\n")
  }
  return(Imat)
}

# initialization

y <- dat[, 1]
X.o <- cbind(1, as.matrix(dat[, -1]))
p <- ncol(X.o)
psi <- 1
M <- 10000

Imat <- SVSS(X.o, y, M)
# filter useful variables
vars1 <- which(colSums(Imat) / M > 0.6) 
# plot(1:M, cumsum(t(Imat[, 7])) / 1:M, type = "l", ylim = c(0, 1))
X.o2 <- model.matrix(y ~ .^2, dat[, vars1])
Imat2 <- SVSS(X.o2, y, M)
vars2 <- colnames(X.o2)[which(colSums(Imat2) / M > 0.9) ]
for (i in vars2) {
  cat(i, ", ", sep = "")
}


pdf("in.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(1:(M-100), cumsum(t(Imat[101:M, 1])) / 1:(M-100), 
     type = "l", ylim = c(0, 1),
     xlab = "Intercept", ylab = "Fraction")
plot(1:(M-100), cumsum(t(Imat[101:M, 4])) / 1:(M-100), 
     type = "l", ylim = c(0, 1),
     xlab = "X3", ylab = "Fraction")
dev.off()

pdf("out.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(1:(M-100), cumsum(t(Imat[101:M, 20])) / 1:(M-100), 
     type = "l", ylim = c(0, 1),
     xlab = "X19", ylab = "Fraction")
plot(1:(M-100), cumsum(t(Imat[101:M, 50])) / 1:(M-100), 
     type = "l", ylim = c(0, 1),
     xlab = "X49", ylab = "Fraction")
dev.off()
