data {
 int<lower=0> n1;
 int<lower=0> n2;
 int y1[n1];
 int y2[n2];
 
} parameters {

  real<lower=0> mu1;
  real<lower=0> mu2;
  real<lower=0> phi1;
  real<lower=0> phi2;
  
} model {
  y1 ~ neg_binomial_2(mu1, phi1);
  y2 ~ neg_binomial_2(mu2, phi2);

  phi1 ~ gamma(0.01, 0.01);
  phi2 ~ gamma(0.01, 0.01);  

} generated quantities{

  int y1rep[n1];
  int y2rep[n2];
  real log_lik[n1+n2];
  for(i in 1:n1){
    y1rep[i] = neg_binomial_2_rng(mu1, phi1);
    log_lik[i] = neg_binomial_2_lpmf(y1[i] | mu1, phi1);
  }
  for(i in 1:n2){
    y2rep[i] = neg_binomial_2_rng(mu2, phi2);
    log_lik[n1+i] = neg_binomial_2_lpmf(y2[i] | mu2, phi2);
  }
}
