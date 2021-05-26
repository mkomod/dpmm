# Dirichlet Process partial linear regression
Rcpp::sourceCpp("./plr.cpp")

# test data
set.seed(1)
n <- 200
K <- 2
x <- seq(0, 4, length.out=n)
b.true <- -2
g.true <- sample(1:K, n, prob=c(0.5, 0.5), replace=T)
eta.true <- c(2, 10)
y <- eta.true[g.true] + x*b + rnorm(n)


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

# sampler
N <- 3e3                               # MCMC number of iters
max.groups <- 10
a.0 <- 1e-60
H0 <- function() rnorm(1, 0, 5)
eta.kern.sd <- 0.3
beta.kern.sd <- 0.4
sigma.kern.sd <- 0.4

# init
g <- matrix(rep(1, n), ncol=1)         # groups are 0 indexed
eta <- matrix(c(rnorm(1), rep(0, max.groups - 1)), ncol=1, nrow=max.groups)
beta <- rnorm(1)
sigma <- rgamma(1, 1, 1)
# sigma <- 1

G <- matrix(0, ncol=N, nrow=n)
ETA <- matrix(0, ncol=N, nrow=max.groups)
BETA <- matrix(0, ncol=N, nrow=1)
SIG <- matrix(0, ncol=N, nrow=1)

# sampler
for (iter in 1:N) {
    # update c
    u <- update_g(g, eta, y, x, beta, sigma, a.0, H0)
    g <- u$g; eta <- u$eta
    
    eta <- update_eta(eta, y, x, g, beta, sigma, eta.kern.sd)
    beta <- update_beta(beta, y, x, eta, g, sigma, beta.kern.sd)
    sigma <- update_sigma(sigma, y, x, eta, g, beta, sigma.kern.sd)

    # save
    G[ , iter] <- g
    ETA[ , iter] <- eta
    BETA[ , iter] <- beta
    SIG[ , iter] <- sigma
}

plot(BETA[1, ])
plot(SIG[1, ])
plot(ETA[2, ])
plot(ETA[3, ])


f1 <- MASS::kde2d(SIG[1, ], BETA[1, ], n=500)
pdf(file="dp_plr_mu_sig.pdf", width=6, height=6)
par(mfrow=c(2,2))
par(mar=c(0,3,3,0))
plot(density(BETA), main=expression(mu), xlab="mu", ylab="", yaxt="n")
par(mar=c(0,0,3,3))
image(f1, main=expression(sigma), col=hcl.colors(12, "purples", rev=T), 
      useRaster=T, xaxt="n", las=1)
par(mar=c(3,3,0,0))
plot.new()
par(mar=c(3,0,0,3))
plot(density(SIG), main="", yaxt="n")
dev.off()
