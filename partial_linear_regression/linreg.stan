// Bayesian Linear Regression with Stan
data {
  int<lower=0> N;
  int<lower=1> K;
  int g[N];
  vector[N] x;
  vector[N] y;
}
parameters {
  vector[K] eta;
  real beta;
  real<lower=0> sigma;
}
model {
  for (i in 1:N) {
    y[i] ~ normal(eta[g[i]] + beta * x[i], sigma);
  }
}
