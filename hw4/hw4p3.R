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

calculate.post <- function(theta0, theta1) {
  
  ## Initialize posterior probabilities
  post <- rep(NA, 1000)
  
  for (i in 1:1000) {
    ss1 <- sum((x.random[1:i] - theta0)^2)
    ss2 <- sum((x.random[1:i] - theta1[1])^2)
    ss3 <- sum((x.random[1:i] - theta1[2])^2)
    ss4 <- sum((x.random[1:i] - theta1[3])^2)
    numerator <- exp(-ss1 / (2 * sig^2))
    denominator <- exp(-ss2 / (2 * sig^2)) / 3 + 
      exp(-ss3 / (2 * sig^2)) / 3 + 
      exp(-ss4 / (2 * sig^2)) / 3
    BF <- numerator / denominator
    
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
par(mar = c(4.1, 5.1, 2.1, 1.1))
plot(post1, xlab = "Iteration", ylab = "1st Eigenvalue",
     main = "", cex.lab = 2, cex.axis = 1.5)
dev.off()