source("dpmm.R")
Rcpp::sourceCpp("./dpmm.cpp")
set.seed(1)

# Test Data, a mixture of three normals
n <- 150
mus <- as.matrix(seq(-3, 3, by=3)) # -3, 0, 3 are the true underlying means
# Draw the group index, note that +1 as cpp index is 0-based and R is 1-based
c.t <- rdiscrete(n, rep(1/3, 3)) + 1
y <- rnorm(n, mus[c.t], 0.3) # 0.3 is the underlying standard deviation

plot(density(y)) # check density looks correct

N <- 5e3
result <- dpmm_gibbs(data = y, N = N) # using mainly default options from function

PHI <- result$PHI
G <- result$G

# results take 1e3 - 1e4 (burn in of 1e3)
grps <- as.numeric(names(table(G[, 5000])))
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
