# Dirichlet Process Mixture Models
# See: Algorithm 7, Neal (2000)

Rcpp::sourceCpp("./dpmm.cpp")
set.seed(1)


# Test Data ---
n <- 150
mus <- as.matrix(seq(-3, 3, by=3))
c.t <- rdiscrete(n, rep(1/3, 3)) + 1   # +1 as cpp index is 0-based and R is 1-based
y <- rnorm(n, mus[c.t], 0.3)

plot(density(y))                       # check density is mixture of normals

# Gibbs sampler w/ MH updates for phi ---
set.seed(1)
g <- matrix(1, n)
phi <- matrix(rnorm(1), length(unique(g)))
g0 <- function() rnorm(1, 0, 3)
a.0 <- 1e-30
kern.sd <- 0.1

N <- 5e3
PHI <- matrix(0, nrow=500, ncol=N)
G <- matrix(0, nrow=n, ncol=N)

for (iter in 1:N) {
    
    # update c (create new groups)
    for (i in 1:n) {
	g.old <- g[i]
	if (any(g[-i] == g[i])) {
	    g.new <- max(g) + 1               # create a new group
	    phi.new <- g0()
	    a <- min(1, a.0/(n-1)*dnorm(y[i], phi.new, 0.3)/dnorm(y[i], phi[g.old], 0.3))
	    if (a > runif(1)) {
		g[i] <- g.new
		phi <- as.matrix(rbind(phi, phi.new))
		colnames(phi) <- NULL
	    }
	} else {
	    tab.g <- table(g[-i])
	    g.new <- sample(as.numeric(names(tab.g)), 1, prob=as.numeric(tab.g))
	    a <- min(1, (n-1)/a.0*dnorm(y[i], phi[g.new],0.3)/dnorm(y[i], phi[g.old], 0.3))
	    g[i] <- ifelse(a > runif(1), g.new, g.old)
	}
    }
    
    # update c
    for (i in 1:n) {
	g.old <- g[i]
	if (any(g[-i] == g[i])) {
	    tab.g <- table(g[-i])
	    g[i] <- sample(unique(g), 1, 
		prob=as.numeric(tab.g)*dnorm(y[i], phi[as.numeric(names(tab.g))], 0.3))
	}
    }
   
    # update phi
    for (i in 1:n) {
	phi.old <- phi[g[i]]
	phi.new <- rnorm(1, phi.old, kern.sd)
	
	p0 <- dnorm(y[i], phi.new, 0.3, log=TRUE) + dnorm(phi.new, 0, 5, log=TRUE) + 
	    dnorm(phi.old, phi.new, kern.sd, log=TRUE) 
	p1 <- dnorm(y[i], phi.old, 0.3, log=TRUE) + dnorm(phi.old, 0, 5, log=TRUE) + 
	    dnorm(phi.new, phi.old, kern.sd, log=TRUE)

	a <- min(1, exp(p0 - p1))

	phi[g[i]] <- ifelse(a > runif(1), phi.new, phi.old)
    }
    
    # save
    PHI[ , iter] <- c(phi, rep(0, 500-length(phi)))
    G[ , iter] <- g
}

# results take 1e3 - 1e4 (burn in of 1e3)
grps <- as.numeric(names(table(g)))
par(mfrow=c(3,3), las=1, mar=c(0,0,0,0))
for(i in seq_along(grps)) {
    gi <- grps[i]
    for (j in seq_along(grps)) {
	gj <- grps[j]
	if (i == j) {
	    plot(density(PHI[gi , 1.5e3:N]), xlab="", ylab="", main="", yaxt="n")
	} else if (i < j) {
	    if (j == i + 1) {
		plot(PHI[gj, 1.5e3:N], PHI[gi, 1.5e3:N], xlab="", ylab="", 
		     col=rgb(0, 0, 1, 0.05), pch=2,
		     xaxt="n", las=1)
	    } else if (j == i + 2) {
		plot(PHI[gj, 1.5e3:N], PHI[gi, 1.5e3:N], xlab="", ylab="", 
		     col=rgb(0, 0, 1, 0.05), pch=2,
		     yaxt="n", las=1)
	    } else {
		plot(PHI[gj, 1.5e3:N], PHI[gi, 1.5e3:N], xlab="", ylab="", 
		     col=rgb(0, 0, 1, 0.05), pch=2,
		     yaxt="n", xaxt="n")
	    }
	} else {
	    plot.new()
	}
    }
}

