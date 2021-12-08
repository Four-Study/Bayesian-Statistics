######### Problem 2 #########

## Import the data and modify column names
dat <- read.csv("Model1_5444.txt", skip = 2, header = FALSE)
colnames(dat) <- c("Y", paste0(rep("X", 50), 1:50))

## augment the data set with interactions
all.x <- model.matrix(Y ~ .^2, dat)
dat.aug <- data.frame(cbind(Y = dat[, 1], all.x[, -1]))

## LASSO
library(glmnet)
lambda <- exp(seq(-10, 1, 0.1))
la.fit1 <- glmnet(dat.aug[, -1], dat.aug[, 1], alpha = 1, lambda = lambda)


plot(la.fit1$lambda, la.fit1$dev.ratio, type = 'l')
plot(la.fit1)
rownames(coef(la.fit1, s = 0.45))[as.numeric(coef(la.fit1, s = 0.45) != 0) == 1]
lm.fit1 <- lm(Y ~ X3 + X17 + X27 + X45 + 
                X3.X15 + X31.X45, data = dat.aug)
summary(lm.fit1)

## AIC
lm.fit2 <- lm(Y ~ X3 + X17 + X27 + X45 + X31.X45, data = dat.aug)
AIC(lm.fit1)
AIC(lm.fit2)

## BIC
BIC(lm.fit1)
BIC(lm.fit2)

## SSVS
