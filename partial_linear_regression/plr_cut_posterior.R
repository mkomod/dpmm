library("rstan")

# test data
set.seed(1)
n <- 200
K <- 2
b <- -2
g <- sample(1:K, n, prob=c(0.5, 0.5), replace=T)
x <- seq(0, 4, length.out=n)
eta.true <- c(2, 10)
y <- eta.true[g] + x*b + rnorm(n)


