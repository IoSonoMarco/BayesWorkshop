library(rstan)
library(bayesplot)

# make data
Y <- c(1, 3)
N <- c(29, 7) # plug 10 instead of 29

### Compile model
model <- stan_model('models/model_binomial.stan')

### Prepare a list of the data
data <- list(y1=Y[1], y2=Y[2], N1=N[1], N2=N[2])

### Sample
fit = sampling(model, data=data, iter=2000, chains=2, warmup=500)

# plots
draws <- extract(fit)
theta1 <- draws$theta1
theta2 <- draws$theta2

plot(density(theta1))
abline(v=data$y1/data$N1, col='red')

plot(density(theta2))
abline(v=data$y2/data$N2, col='red')

# Inspect bivariate posteriors
mcmc_pairs(fit, pars=c("theta1", "theta2"), diag_fun="dens")
mcmc_areas(fit, pars=c('theta1', 'theta2', 'theta_diff'), prob = 0.95)
