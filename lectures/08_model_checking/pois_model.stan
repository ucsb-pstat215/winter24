data {
 int<lower=0> n1;
 int<lower=0> n2;
 int y1[n1];
 int y2[n2];
 
} parameters {

  real<lower=0> lambda1;
  real<lower=0> lambda2;
  
} model {
  y1 ~ poisson(lambda1);
  y2 ~ poisson(lambda2);

} generated quantities{

  int y1rep[n1];
  int y2rep[n2];
  real log_lik[n1+n2];
  for(i in 1:n1){
    y1rep[i] = poisson_rng(lambda1);
    log_lik[i] = poisson_lpmf(y1[i] | lambda1);
  }
  for(i in 1:n2){
    y2rep[i] = poisson_rng(lambda2);
    log_lik[n1+i] = poisson_lpmf(y2[i] | lambda2);
  }
}
