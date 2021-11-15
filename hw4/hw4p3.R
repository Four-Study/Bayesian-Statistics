## hyper-parameter
sig <- sqrt(5)
set.seed(12345)
x <- c(rnorm(400, 10, sig), 
       rnorm(200, 15, sig), 
       rnorm(300, 17, sig), 
       rnorm(100, 20, sig))
hist(x, breaks = 20)
