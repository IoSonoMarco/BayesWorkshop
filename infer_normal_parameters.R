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
