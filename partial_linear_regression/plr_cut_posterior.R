Rcpp::sourceCpp("./plr_cut.cpp")

# test data
set.seed(1)
n <- 200
K <- 2
x <- seq(0, 4, length.out=n)
b.true <- -2
g.true <- sample(1:K, n, prob=c(0.5, 0.5), replace=T)
eta.true <- c(2, 10)
y <- eta.true[g.true] + x*b.true + rnorm(n)


update_eta <- function(eta, y, x, g, beta, sigma, eta.kern.sd) 
{
    for (i in 1:n) {
	eta.old <- eta[g[i]]
	eta.new <- rnorm(1, eta.old, eta.kern.sd)
	
	p0 <- log_marginal_lik(y[i], x[i], eta.new, sigma, sigma_0)
	    dnorm(eta.new, 0, 5, log=TRUE)

	p1 <- log_marginal_lik(y[i], x[i], eta.old, sigma, sigma_0)
	    dnorm(eta.old, 0, 5, log=TRUE)

	a <- min(1, exp(p0 - p1))
	eta[g[i]] <- ifelse(a > runif(1), eta.new, eta.old)
    }
    return(eta)
}

N <- 3e3                               # MCMC iters
max.groups <- 10                       # set an upper limit on the number of groups
G.max <- 1                             # count for groups
a.0 <- 1e-30                           # concentration rate
H0 <- function() rnorm(1, 0, 5)        # base measure for DP
eta.kern.sd <- 0.3                     # kernel sd for eta

# init vars
g <- matrix(rep(1, n), ncol=1)         # init to same group
eta <- matrix(c(rnorm(1), rep(0,max.groups-1)), ncol=1, nrow=max.groups)
sigma_0 <- .1                             # prior variance for beta
sigma <- 1

ETA <- matrix(0, ncol=N, nrow=max.groups)
G <- matrix(0, ncol=N, nrow=n)


for (iter in 1:N) {
    # update c
    u <- update_g_marginal(g, eta, y, x, sigma, sigma_0, a.0, H0)
    g <- u$g; eta <- u$eta
    
    # update eta
    eta <- update_eta(eta, y, x, g, beta, sigma, eta.kern.sd)

    # save
    G[ , iter] <- g
    ETA[ , iter] <- eta
}

table(g)
mean(ETA[1, ]); sd(ETA[1, ])
mean(ETA[2, ]); sd(ETA[2, ])

cov(ETA[1, ], ETA[2,])
mean(ETA[1, ]) - mean(ETA[2, ])

pdf("marginal_data.pdf", width=6, height=5)
plot(c(0, 4), c(4, 4), type="l", lwd=2, col="orange", ylim=c(-7, 8), main="", ylab="y", xlab="x" )
lines(c(0, 4), c(-4, -4), type="l", lwd=2, col="darkorchid1")
y.centered <- eta.true[g.true] -6 + rnorm(n)
cls <- c("darkorchid1", "orange")
points(x, y.centered, pch=20, col=cls[g.true])
dev.off()


s1 <- stan("./linreg.stan", data=list(N=n, K=2, g=as.numeric(g), x=x, y=y))
samples <- extract(s1, c("eta[1]", "eta[2]", "beta", "sigma"))
X <- rbind(samples$`eta[1]`, samples$`eta[2]`, samples$beta, samples$sigma)

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

