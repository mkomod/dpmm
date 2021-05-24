# Dirichlet Process Mixture Models
# See: Algorithm 7, Neal (2000)

Rcpp::sourceCpp("./dpmm.cpp")
set.seed(1)

# Test Data, a mixture of three normals
n <- 150
mus <- as.matrix(seq(-3, 3, by=3)) # -3, 0, 3 are the true underlying means
# Draw the group index, note that +1 as cpp index is 0-based and R is 1-based
c.t <- rdiscrete(n, rep(1/3, 3)) + 1
y <- rnorm(n, mus[c.t], 0.3) # 0.3 is the underlying standard deviation

plot(density(y)) # check density looks correct

# Gibbs sampler w/ MH updates for phi
set.seed(1)
g <- matrix(1, n) # init group membership, all in first group
phi <- matrix(rnorm(1), length(unique(g))) # init phi randomly
g0 <- function() rnorm(1, 0, 3)
a.0 <- 1e-30
kern.sd <- 0.1

N <- 5e3 # the number of MCMC iterations
G.max <- 1 # init the largest numbered group at 1
PHI <- matrix(0, nrow=500, ncol=N) # assuming maximum number of groups <500
G <- matrix(0, nrow=n, ncol=N)

for (iter in 1:N) {
    # Update c (create new groups)
    for (i in 1:n) {
	    g.old <- g[i]
	    # When c[i] is not a singleton, such that c[i] = c[j] for some j
	    if (any(g[-i] == g[i])) {
	      g.new <- G.max + 1 # create a new group
	      phi.new <- g0() # draw parameters from g0
	      r <- a.0/(n-1)*dnorm(y[i], phi.new, 0.3)/dnorm(y[i], phi[g.old], 0.3)
	      a <- min(1, r)
	      if (a > runif(1)) {
		      g[i] <- g.new
		      G.max <- G.max + 1 # update G.max
		      phi <- as.matrix(rbind(phi, phi.new))
		      rownames(phi) <- NULL
	      }
	  } else {
	    # When c[i] is a singleton
	    tab.g <- table(g[-i])
	    # Sample new group proportional to existing groups
	    g.new <- sample(as.numeric(names(tab.g)), 1, prob=as.numeric(tab.g))
	    r <- (n-1)/a.0*dnorm(y[i], phi[g.new],0.3)/dnorm(y[i], phi[g.old], 0.3)
	    a <- min(1, r)
	    g[i] <- ifelse(a > runif(1), g.new, g.old)
	  }
  }
    
  # Update c
  for (i in 1:n) {
    # Only for non-singletons
	  if (any(g[-i] == g[i])) {
	    tab.g <- table(g[-i])
	    g.vals <- as.numeric(names(tab.g))
	    g[i] <- sample(g.vals, 1, prob=as.numeric(tab.g)*dnorm(y[i], phi[g.vals], 0.3))
	  }
  }
   
  # Update phi
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
    
  # Save
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
