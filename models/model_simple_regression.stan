data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}

parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma2;
}

transformed parameters {
  real<lower=0> sigma;
  sigma = sqrt(sigma2);
}

model {
  sigma2 ~ inv_gamma(1, 1);
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  for(n in 1:N) y[n] ~ normal(beta0 + beta1*x[n], sigma);
}

generated quantities {
  vector[N] log_lik; 
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | beta0 + beta1*x[n], sigma);
}
