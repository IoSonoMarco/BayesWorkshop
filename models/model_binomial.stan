data {
  int<lower=1> N1;
  int<lower=1> N2;
  int<lower=1> y1;
  int<lower=1> y2;
}

parameters {
  real<lower=0, upper=1> theta1;
  real<lower=0, upper=1> theta2;
}

transformed parameters {
  real theta_diff;
  theta_diff = theta1 - theta2;
}

model {
  y1 ~ binomial(N1, theta1);
  y2 ~ binomial(N2, theta2);
  theta1 ~ beta(1,1);
  theta2 ~ beta(1,1);
}

