data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] x;
  int<lower=0> I[N];
}

parameters {
  real mu;
  real<lower=0> tau2;
  real<lower=0> sigma2;
  vector[K] beta;
} 

transformed parameters {
  real<lower=0> tau;
  real<lower=0> sigma;
  tau = sqrt(tau2);
  sigma = sqrt(sigma2);
}

model {

  for (n in 1:N) { 
    x[n] ~ normal(beta[I[n]], sigma);
  }
  
  mu ~ normal(0, 10.0);
  tau2 ~ inv_gamma(0.5, 0.5);
  sigma2 ~ inv_gamma(0.5, 0.5);
  beta ~ normal(mu, tau);
  

  
}
