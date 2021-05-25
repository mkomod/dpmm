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


likelihood <- function(y.i, x, eta, beta, sigma) {
    dnorm(y.i, eta + beta * x, sigma)
}

log_likelihood <- function(y, x, eta, beta, sigma) {
    sum(dnorm(y, eta + beta*x, sigma, log=T))
}


# sampler
N <- 5e3                               # MCMC iters
max.groups <- n                        # set an upper limit on the number of groups
G.max <- 1                             # count for groups
a.0 <- 1e-20                           # concentration rate
H0 <- function() rnorm(1, 0, 5)        # base measure for DP
eta.kern.sd <- 0.3                     # kernel sd for eta

# init vars
g <- matrix(rep(1, n), ncol=1)         # init to same group
eta <- matrix(c(rnorm(1), rep(0, 24)), ncol=1, nrow=max.groups)
beta <- -2                             # equiv to dirac mass on beta
sigma <- 1                             # equiv to dirac mass on sigma

ETA <- matrix(0, ncol=N, nrow=n)
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
	if (any(g[-i] == g[i])) {
	    gs <- table(g[-i])
	    g.i <- as.numeric(names(gs))
	    g.p <- as.numeric(gs)
	    if (length(g.p) == 1) {
		g[i] <- g.i  		# if there's only one group set it to that
	    } else {
		g[i] <- sample(g.i, 1, 
		    prob=g.p*likelihood(y[i], x[i], eta[g.i], beta, sigma))
	    }
	}
    }

    # update eta - MH within Gibbs
    for (i in 1:n) {
	eta.old <- eta[g[i]]
	eta.new <- rnorm(1, eta.old, eta.kern.sd) # note: our kernel is symmetric
	
	p0 <- log_likelihood(y[i], x[i], eta.new, beta, sigma) +
	    dnorm(eta.new, 0, 5, log=TRUE)

	p1 <- log_likelihood(y[i], x[i], eta.old, beta, sigma) +
	    dnorm(eta.old, 0, 5, log=TRUE)

	a <- min(1, exp(p0 - p1))
	eta[g[i]] <- ifelse(a > runif(1), eta.new, eta.old)
    }

    # save
    G[ , iter] <- g
    ETA[ , iter] <- eta
}


mean(ETA[1, ])
mean(ETA[2, ])

f1 <- MASS::kde2d(ETA[1, 1e3:5e3], ETA[2, 1e3:5e3], n=500)
pdf("eta_post.pdf", width=6, height=5)
image(f1, col=hcl.colors(12), xlab=expression(eta[1]), ylab=expression(eta[2]),
las=1, useRaster=TRUE)
dev.off()


cls <- c("darkorchid1", "orange")
pdf(file="plr.pdf", width=6, height=5)
plot(x, y, pch=20, col=cls[g], las=1)
dev.off()

