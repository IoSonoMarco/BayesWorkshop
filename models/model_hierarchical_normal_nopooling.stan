data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] x;
  int<lower=0> I[N];
}

parameters {
  real<lower=0> sigma2;
  vector[K] beta;
} 

transformed parameters {
  real<lower=0> sigma;
  sigma = sqrt(sigma2);
}

model {
    
  for (n in 1:N) {
      x[n] ~ normal(beta[I[n]], sigma);
  }

  sigma2 ~ inv_gamma(0.5, 0.5);
  beta ~ normal(0, 10.0);
  
}
