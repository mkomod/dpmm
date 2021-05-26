#' MCMC sampler for Dirichlet Process Mixture Models (DPMM)
#' 
#' Implements the Gibbs sampler (with Metropolis-Hastings step for phi) for
#' DPMMs, following Algorithm 7, Neal (2000).
#'
#' @param data A numeric vector.
#' @param g0
#' @param a.0
#' @param kern.sd
#' @param N The number of MCMC iterations.
dpmm_gibbs <- function(data, 
                       g0 = function() rnorm(1, 0, 3), 
                       a.0 = 1e-30, 
                       kern.sd = 0.1,
                       N = 5e3) {
  
  set.seed(1)
  n <- length(data)
  
  # Initialisation and MCMC storage
  g <- matrix(1, n) # group membership, all in first group
  phi <- matrix(rnorm(1), length(unique(g))) # phi, at random
  G.max <- 1 # the largest numbered group
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
  
  return(list(PHI = PHI, G = G))
}
