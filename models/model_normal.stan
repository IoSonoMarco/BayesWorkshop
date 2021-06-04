data {
  int<lower=1> N;
  vector[N] x;
}

parameters {
  real mu;
  real<lower=0> sigma2;
}

transformed parameters {
  real<lower=0> sigma;
  sigma = sqrt(sigma2);
}

model {
  mu ~ normal(0, 10);
  sigma2 ~ inv_gamma(1, 1);
  x ~ normal(mu, sigma);
}

