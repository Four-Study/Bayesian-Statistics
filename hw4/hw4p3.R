## hyper-parameter
sig <- sqrt(5)
set.seed(12345)

## generate the data
x <- c(rnorm(400, 10, sig), 
       rnorm(200, 15, sig), 
       rnorm(300, 17, sig), 
       rnorm(100, 20, sig))

## draw a histogram for x
pdf("imgs/hist_x.pdf", width = 10, height = 10)
par(mar = c(4.1, 5.1, 2.1, 1.1))
hist(x, breaks = 20, main = "", 
     cex.lab = 2, cex.axis = 1.5)
dev.off()

## shuffle x and draw a scatter plot
set.seed(202111)
x.random <- sample(x, length(x))

pdf("imgs/scatter_x.pdf", width = 10, height = 10)
par(mar = c(4.1, 5.1, 2.1, 1.1))
plot(x.random, main = "", ylab = "Shuffled x",
     cex.lab = 2, cex.axis = 1.5)
dev.off()

## write a function to calculate the posterior probability

library(matrixStats)
calculate.post <- function(theta0, theta1) {
  
  ## Initialize posterior probabilities
  post <- rep(NA, 1000)
  
  for (i in 1:1000) {
    ss1 <- sum((x.random[1:i] - theta0)^2) / (2 * sig^2)
    ss2 <- sum((x.random[1:i] - theta1[1])^2) / (2 * sig^2)
    ss3 <- sum((x.random[1:i] - theta1[2])^2) / (2 * sig^2)
    ss4 <- sum((x.random[1:i] - theta1[3])^2) / (2 * sig^2)
    logBF <- -ss1 - logSumExp(-c(ss2, ss3, ss4))
    BF <- 3 * exp(logBF)
    
    post[i] <- 1 / (1 + 3 * 1 / BF)
  }
  return(post)
}

post1 <- calculate.post(10, c(15, 17, 20))
post2 <- calculate.post(15, c(10, 17, 20))
post3 <- calculate.post(17, c(10, 15, 20))
post4 <- calculate.post(20, c(10, 15, 17))

## draw 
pdf("imgs/posterior.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(4.1, 2.5, 2.1, 1.1))
plot(post1, xlab = "Number of Points", ylab = "",
     main = "H0 = 10", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
     type = "l", ylim = c(0, 1))
plot(post2, xlab = "Number of Points", ylab = "",
     main = "H0 = 15", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
     type = "l", ylim = c(0, 1))
plot(post3, xlab = "Number of Points", ylab = "",
     main = "H0 = 17", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
     type = "l", ylim = c(0, 1))
plot(post4, xlab = "Number of Points", ylab = "",
     main = "H0 = 20", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
     type = "l", ylim = c(0, 1))
dev.off()

## write a function to calculate p-value in sequence

calculate.pvalue <- function(theta0) {
  
  ## Initialize posterior probabilities
  pvalue <- rep(NA, 1000)
  
  for (i in 1:1000) {
    xbar <- mean(x.random[1:i])
    pvalue[i] <- 2 * min(pnorm(xbar, theta0, sig / sqrt(i)), 
                         1-pnorm(xbar, theta0, sig / sqrt(i)))
  }
  return(pvalue)
}

pvalue1 <- calculate.pvalue(10)
pvalue2 <- calculate.pvalue(15)
pvalue3 <- calculate.pvalue(17)
pvalue4 <- calculate.pvalue(20)

pdf("imgs/pvalue.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(4.1, 2.5, 4.1, 1.1))
plot(pvalue1, xlab = "Number of Points", ylab = "",
     main = "H0 = 10", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
     type = "l", ylim = c(0, 1))
plot(pvalue2, xlab = "Number of Points", ylab = "",
     main = "H0 = 15", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
     type = "l", ylim = c(0, 1))
plot(pvalue3, xlab = "Number of Points", ylab = "",
     main = "H0 = 17", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
     type = "l", ylim = c(0, 1))
plot(pvalue4, xlab = "Number of Points", ylab = "",
     main = "H0 = 20", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,
     type = "l", ylim = c(0, 1))
dev.off()

save.image("p3.RData")
