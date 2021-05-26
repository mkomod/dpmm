update_eta <- function(eta, y, x, g, beta, sigma, eta.kern.sd) 
{
    for (i in 1:n) {
	eta.old <- eta[g[i]]
	eta.new <- rnorm(1, eta.old, eta.kern.sd)
	
	p0 <- dp_plr_log_lik(y[i], x[i], eta.new, beta, sigma) +
	    dnorm(eta.new, 0, 5, log=TRUE)

	p1 <- dp_plr_log_lik(y[i], x[i], eta.old, beta, sigma) +
	    dnorm(eta.old, 0, 5, log=TRUE)

	a <- min(1, exp(p0 - p1))
	eta[g[i]] <- ifelse(a > runif(1), eta.new, eta.old)
    }
    return(eta)
}


update_beta <- function(beta, y, x, eta, g, sigma, beta.kern.sd) 
{
    for (i in 1:n) {
	beta.old <- beta
	beta.new <- rnorm(1, beta.old, beta.kern.sd)
	
	p0 <- dp_plr_log_lik(y[i], x[i], eta[g[i]], beta.new, sigma) +
	    dnorm(beta.new, 0, 5, log=TRUE)

	p1 <- dp_plr_log_lik(y[i], x[i], eta[g[i]], beta.old, sigma) +
	    dnorm(beta.old, 0, 5, log=TRUE)

	a <- min(1, exp(p0 - p1))
	beta <- ifelse(a > runif(1), beta.new, beta.old)
    }
    return(beta)
}


update_sigma <- function(sigma, y, x, eta, g, beta, sigma.kern.sd) 
{
    for (i in 1:n) {
	sigma.old <- sigma
	sigma.new <- rlnorm(1, sigma.old, sigma.kern.sd)
	
	p0 <- dp_plr_log_lik(y[i], x[i], eta[g[i]], beta, sigma) +
	    dlnorm(sigma.old, sigma.new, sigma.kern.sd, log=T)

	p1 <- dp_plr_log_lik(y[i], x[i], eta[g[i]], beta, sigma) +
	    dlnorm(sigma.new, sigma.old, sigma.kern.sd, log=T)

	a <- min(1, exp(p0 - p1))
	sigma <- ifelse(a > runif(1), sigma.new, sigma.old)
    }
    return(sigma)
}

g.old <- g
N <- 3e3                               # MCMC number of iters
eta.kern.sd <- 0.3
beta.kern.sd <- 0.4
sigma.kern.sd <- 0.4

# init
beta <- rnorm(1)
sigma <- rgamma(1, 1, 1)

ETA <- matrix(0, ncol=N, nrow=max.groups)
BETA <- matrix(0, ncol=N, nrow=1)
SIG <- matrix(0, ncol=N, nrow=1)

for (iter in 1:N) {
    # update c
    eta <- update_eta(eta, y, x, g, beta, sigma, eta.kern.sd)
    beta <- update_beta(beta, y, x, eta, g, sigma, beta.kern.sd)
    sigma <- update_sigma(sigma, y, x, eta, g, beta, sigma.kern.sd)

    # save
    ETA[ , iter] <- eta
    BETA[ , iter] <- beta
    SIG[ , iter] <- sigma
}


X <- rbind(ETA[c(1,2), ], BETA, SIG)
v <- c(expression(eta[1]), expression(eta[2]), expression(mu), expression(sigma))

pdf("cut_post.pdf", width=8, height=8)
par(mfrow=c(4,4))
for (i in 1:4) {
    for (j in 1:4) {
	if (i == 1 && j == 1) par(mar=c(0, 2, 0, 0))
	if (i == 1 && j > 1 && j < 4) par(mar=c(0, 0, 0, 0))
	if (i > 1 && i < 4 && j > 1 && j <= 4) par(mar=c(0, 0, 0, 0))

	if (i == j) {
	    plot(density(X[i, ]), main="", yaxt="n")
	    mtext(v[i], 2, .5, las=1)
	} else if (j > i) {
	    f <- MASS::kde2d(X[j, ], X[i, ], n=200)
	    if(j == 4 && i == 3) {
	    image(f, useRaster=T, las=1,  col=hcl.colors(12, "purples", rev=T))
	    } else if (j == i + 1) {
	    image(f, useRaster=T, las=1,  col=hcl.colors(12, "purples", rev=T), xaxt="n")
	    } else {
	    image(f, useRaster=T, las=1,  col=hcl.colors(12, "purples", rev=T), xaxt="n", yaxt="n")
	    }
	} else {
	    plot.new()
	}
    }
}
dev.off()
