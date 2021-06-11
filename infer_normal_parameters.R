library(rstan)
library(bayesplot)

### Simulate data
N <- 50
mu <- 100
sigma <- 15
x <- rnorm(N, mu, sigma)

### Compile model
model <- stan_model('models/model_normal.stan')

### Prepare a list of the data
data <- list(x=x, N=N)

### Sample
fit = sampling(model, data=data, iter=2000, chains=2, warmup=500)

draws <- extract(fit)
sigma2 <- draws$sigma2
sigma <- draws$sigma

plot(density(sigma))
plot(density(sqrt(sigma2)))


# Inspect summaries and R_hat diagnostics
print(fit)


# Inspect bivariate posteriors
mcmc_pairs(fit, pars=c("mu", "sigma"), diag_fun="dens")

mu <- draws$mu

draws <- extract(fit)
sim <- draws$sim

hist(x)
hist(sim)

hist(x, freq=F)
lines(density(sim))

## sim ##

new_x <- seq(-150,150,length=1000)
hist(x, freq=F)
for(iter in 1:100){
  sample_pos <- sample(3000, 1)
  mu_sample <- mu[sample_pos]
  sigma_sample <- sigma[sample_pos]
  lines(new_x, dnorm(new_x, mu_sample, sigma_sample), type='l', col='blue')
}
