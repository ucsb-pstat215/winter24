data {
 int<lower=0> n1;
 int<lower=0> n2;
 int y1[n1];
 int y2[n2];
 
} parameters {

  real<lower=0> lambda1;
  real<lower=0> lambda2;
  real<lower=0, upper=1> theta1;
  real<lower=0, upper=1> theta2;
  
} model {

  for (n in 1:n1) {
    if (y1[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | theta1),
                            bernoulli_lpmf(0 | theta1)
                              + poisson_lpmf(y1[n] | lambda1));
    else
      target += bernoulli_lpmf(0 | theta1)
                  + poisson_lpmf(y1[n] | lambda1);
  }
  for (n in 1:n2) {
    if (y2[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | theta2),
                            bernoulli_lpmf(0 | theta2)
                              + poisson_lpmf(y2[n] | lambda2));
    else
      target += bernoulli_lpmf(0 | theta2)
                  + poisson_lpmf(y2[n] | lambda2);
  }
} generated quantities{

  int y1rep[n1];
  int y2rep[n2];
  real log_lik[n1+n2];
  
  for(i in 1:n1){
    if(bernoulli_rng(theta1)==1) {
      y1rep[i] = 0;
    }
    else {
      y1rep[i] = poisson_rng(lambda1);
    }
    if(y1[i] == 0)
      log_lik[i] = log(theta1 + (1-theta1)*exp(-lambda1));
    else
      log_lik[i] = log((1-theta1)) + poisson_lpmf(y1[i] | lambda1);
  }
  for(i in 1:n2){
    if(bernoulli_rng(theta2)==1){
      y2rep[i] = 0;
    }
    else {
      y2rep[i] = poisson_rng(lambda2);
    }
    if(y2[i] == 0)
      log_lik[i+n1] = log(theta2 + (1-theta2)*exp(-lambda2));
    else
      log_lik[i+n1] = log((1-theta2)) + poisson_lpmf(y2[i] | lambda2);
  }
}