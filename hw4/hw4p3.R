## hyper-parameter
sig <- sqrt(5)
set.seed(12345)

## generate the data
x <- c(rnorm(400, 10, sig), 
       rnorm(200, 15, sig), 
       rnorm(300, 17, sig), 
       rnorm(100, 20, sig))

## draw a histogram for x
pdf("hist_x.pdf", width = 10, height = 10)
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

calculate.post <- function(H0, H1) {
  ss <- (x.random - H0)^2
  numerator <- exp()
}