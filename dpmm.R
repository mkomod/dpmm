# dirichlet process mixute model
library(Rcpp)
Rcpp::sourceCpp("./dpmm.cpp")

set.seed(1)

# test cpp funcs
rdirchlet(10, 1:5); apply(rdirchlet(10, 1:5), 1, sum)
rdiscrete(5, rep(1,5)/5); table(rdiscrete(1e4, rep(1,5)/5)) / 1e4


# test data
n <- 150
mus <- as.matrix(seq(-3, 3, by=3))
c.t <- rdiscrete(n, rep(1/3, 3)) + 1   # +1 as cpp index is 0-based and R is 1-based
y <- rnorm(n, mus[c.t], 0.3)

plot(density(y))                       # check density is mixture of normals


# Gibbs sampler for DPMM
# See: Algorithm 7, Neal (2000)

# init vars
g <- matrix(1:2, n)
phi <- matrix(rnorm(2), length(unique(g)))
g0 <- function() rnorm(1, 0, 3)
a <- 0.3

PHI <- list()
G <- matrix(0, nrow=n, ncol=1e4)

for (iter in 1:5e3) {

    for (i in 1:n) {
	g.old <- g[i]
	if (any(g[-i] == g[i])) {
	    g.new <- max(g) + 1               # create a new group
	    phi.new <- g0()
	    a <- min(1, a/(n-1)*dnorm(y[i], phi.new, 0.3)/dnorm(y[i], phi[g.old], 0.3))
	    if (a > runif(1)) {
		g[i] <- g.new
		phi <- as.matrix(rbind(phi, phi.new))
	    }
	} else {
	    g.new <- sample(unique(g[-i]), 1, prob=table(g[-i]))
	    a <- min(1, (n-1)/a*dnorm(y[i], phi[g.new], 0.3)/dnorm(y[i], phi[g.old], 0.3))
	    g[i] <- ifelse(a > runif(1), g.new, g.old)
	}
    }

    for (i in 1:n) {
	g.old <- g[i]
	if (any(g[-i] == g[i])) {
	    g[i] <- sample(unique(g), 1, 
			   prob=table(g[-i])*dnorm(y[i], phi[unique(g[-i])], 0.3))
	}
    }

    for (i in seq_along(unique(g))) {
	phi.old <- phi[i]
	phi.new <- rnorm(1, phi.old, 0.01)
	
	p0 <- sum(dnorm(y, phi.new, 0.3, log=TRUE)) + dnorm(phi.new, log=TRUE) +
	    dnorm(phi.old, phi.new, 0.01, log=TRUE) 
	p1 <- sum(dnorm(y, phi.old, 0.3, log=TRUE)) + dnorm(phi.old, log=TRUE) +
	    dnorm(phi.new, phi.old, 0.01, log=TRUE)

	a <- min(1, exp(p0 - p1))

	phi[i] <- ifelse(a > runif(1), phi.new, phi.old)
    }
    
    PHI[[iter]] <- phi
    G[ , iter] <- g
}


