data {
  int<lower=1> K; // number of mixture components
  int<lower=1> N; // number of data points
  real y[N]; // observations
}
parameters {
  simplex[K] theta; // mixing proportions
  real mu[K]; // locations of mixture components
  real<lower=0> sigma[K]; // scales of mixture components
}
model {
  real ps[K]; // temp for log component densities
  sigma ~ cauchy(0, 1);
  mu ~ normal(10, 1);
  for (n in 1:N) {
    for (k in 1:K) {
      ps[k] = log(theta[k])
      + normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    target += log_sum_exp(ps);
  }
} generated quantities{
  real log_lik[N];
  real ps[K];
  for(n in 1:N) {
    log_lik[n] = 0;
    for (k in 1:K) {
      ps[k] = log(theta[k])
      + normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    log_lik[n] += log_sum_exp(ps);
  }
}