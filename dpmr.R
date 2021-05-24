# Dirichlet Process mixutre of regression models
set.seed(1)

# test data
n <- 200
K <- 2
b <- -2
g <- sample(1:K, n, prob=c(0.5, 0.5), replace=T)
x <- seq(0, 4, length.out=n)
eta.true <- c(2, 10)
y <- eta.true[g] + x*b + rnorm(n)

plot(x, y)

likelihood <- function(y.i, x, eta, beta, sigma) {
    dnorm(y.i, eta + beta * x, sigma)
}

log_likelihood <- function(y, x, eta, beta, sigma) {
    sum(dnorm(y, eta + beta*x, sigma, log=T))
}


# sampler
N <- 2e3                               # MCMC number of iters
max.groups <- n
G.max <- 1
a.0 <- 1e-20
H0 <- function() rnorm(1, 0, 5)
eta.kern.sd <- 0.3
beta.kern.sd <- 0.15
sigma.kern.rate <- 8
sigma.kern.scale <- 8

# init
g <- matrix(rep(1, n), ncol=1)
eta <- matrix(c(rnorm(1), rep(0, 24)), ncol=1, nrow=max.groups)
# beta <- rnorm(1)
beta <- -2
# sigma <- matrix(rgamma(1, 1, 1), ncol=1)
sigma <- 1

ETA <- matrix(0, ncol=N, nrow=n)
BETA <- matrix(0, ncol=N, nrow=1)
SIGMA <- matrix(0, ncol=N, nrow=1)
G <- matrix(0, ncol=N, nrow=n)


for (iter in 1:N) 
{
    # update c
    for (i in 1:n) {
	if (G.max < max.groups && any(g[i] == g[-i])) {
	    g.new <- G.max + 1
	    eta.new <- H0()
	    a <- min(1, a.0 /(n - 1) * 
		     likelihood(y[i], x[i], eta.new,   beta, sigma) /
		     likelihood(y[i], x[i], eta[g[i]], beta, sigma))
	    if (a > runif(1)) {
		g[i] <- g.new
		eta[g[i]] <- eta.new
		G.max <- G.max + 1
	    }
	} else {
	    gs <- table(g[-i])
	    g.i <- as.numeric(names(gs))
	    g.p <- as.numeric(gs)
	    if (length(g.i) == 1) {
		g[i] <- g.i
	    } else {
		g.new <- sample(g.i, 1, prob=g.p)
		a <- min(1, (n - 1)/a.0 * 
			 likelihood(y[i], x[i], eta[g.new], beta, sigma) /
			 likelihood(y[i], x[i], eta[g[i]], beta, sigma))
		g[i] <- ifelse(a > runif(1), g.new, g[i])
	    }
	}
    }
    
    for (i in 1:n) {
	if (any(g[i] == g[-i])) {
	    gs <- table(g[-i])
	    g.i <- as.numeric(names(gs))
	    g.p <- as.numeric(gs)
	    if (length(g.p) == 1) {
		g[i] <- g.i
	    } else {
		g[i] <- sample(g.i, 1, 
		    prob=g.p*likelihood(y[i], x[i], eta[g.i], beta, sigma))
	    }
	}
    }

    # update eta
    for (i in 1:n) {
	eta.old <- eta[g[i]]
	eta.new <- rnorm(1, eta.old, eta.kern.sd)

	p0 <- log_likelihood(y[i], x[i], eta.new, beta, sigma) +
	    dnorm(eta.new, 0, 5, log=TRUE) + 
	    dnorm(eta.old, eta.new, eta.kern.sd, log=TRUE) 

	p1 <- log_likelihood(y[i], x[i], eta.old, beta, sigma) +
	    dnorm(eta.old, 0, 5, log=TRUE) + 
	    dnorm(eta.new, eta.old, eta.kern.sd, log=TRUE) 

	a <- min(1, exp(p0 - p1))

	eta[g[i]] <- ifelse(a > runif(1), eta.new, eta.old)
    }

    # update beta
    # beta.old <- beta
    # beta.new <- rnorm(1, beta.old, beta.kern.sd)
    # p0 <- log_likelihood(y, x, eta[g], beta.new, sigma) +
	    # dnorm(beta.new, 0, 5, log=TRUE) + 
	    # dnorm(beta.old, beta.new, beta.kern.sd, log=TRUE) 

    # p1 <- log_likelihood(y, x, eta[g], beta.old, sigma) +
	    # dnorm(beta.old, 0, 5, log=TRUE) + 
	    # dnorm(beta.new, beta.old, beta.kern.sd, log=TRUE) 

    # a <- min(1, exp(p0 - p1))
    # beta <- ifelse(a > runif(1), beta.new, beta.old)


    # update sigma
    # sigma.old <- sigma
    # sigma.new <- rgamma(1, sigma.old*sigma.kern.scale, rate=sigma.kern.rate)
    # p0 <- log_likelihood(y, x, eta[g], beta, sigma.new) +
	    # dunif(sigma.new, 0, 5, log=TRUE) + 
	    # dgamma(sigma.old, sigma.new*sigma.kern.scale, 
		   # sigma.kern.rate, log=TRUE) 

    # p1 <- log_likelihood(y, x, eta[g], beta, sigma.old) +
	    # dunif(sigma.old, 0, 5, log=TRUE) + 
	    # dgamma(sigma.new, sigma.old*sigma.kern.scale, 
		   # sigma.kern.rate, log=TRUE) 

    # a <- min(1, exp(p0 - p1))
    # sigma <- ifelse(a > runif(1), sigma.new, sigma.old)
    
    # save
    G[ , iter] <- g
    ETA[ , iter] <- eta
    BETA[ , iter] <- beta
    SIGMA[ , iter] <- sigma
}


plot(ETA[1 , ])
plot(ETA[2 , ])
# plot(BETA[1, ])
# plot(SIGMA[1, ])

