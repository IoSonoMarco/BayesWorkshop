library(rstan)
library(loo)
library(bridgesampling)
library(bayesplot)

standardize <- function(x){(x-mean(x))/sd(x)}


data(mtcars)
d <- data.frame(x=mtcars$disp, y=mtcars$mpg)

plot(y ~ x, data=d)

##############
### linear ###
##############

### Stan models
model_linear = stan_model('models/model_simple_regression.stan')

### Stan data
data_linear = list(x = standardize(d$x), 
                   y = standardize(d$y), 
                   N = nrow(d))

data_linear = list(x = d$x, 
                   y = d$y, 
                   N = nrow(d))

### Sample from models
fit_linear = sampling(model_linear, data = data_linear, chains=3, warmup=800, iter=3000)

### Check some plots and predictions
print(fit_linear)

mcmc_pairs(fit_linear, pars=c('beta0','beta1'), diag_fun='dens')

draws <- as.matrix(fit_linear)
mcmc_areas(draws, pars=c('beta1'), prob=0.95)

#################
### quadratic ###
#################

### Stan models
model_quadratic = stan_model('models/model_quadratic_regression.stan')

### Stan data
data_quadratic = list(x = standardize(d$x), 
                      xsq = standardize(d$x^2), 
                      y = standardize(d$y),
                      N = nrow(d))

### Sample from models
fit_quadratic = sampling(model_quadratic, data = data_quadratic, chains=3, warmup=800, iter=3000)

### Check some plots and predictions
print(fit_quadratic)

mcmc_pairs(fit_quadratic, pars=c('beta0','beta1','beta2'), diag_fun='dens')

draws <- as.matrix(fit_quadratic)
mcmc_areas(draws, pars=c('beta1', 'beta2'), prob=0.95)

mcmc_acf_bar(draws, pars=c('beta0','beta1','beta2'))

### Posterior predictive model comparison

# leave-one-out
loo_linear <- loo(fit_linear)
loo_quadratic <- loo(fit_quadratic)

# Inspect comparison
loo_compare(loo_linear, loo_quadratic)

### Model comparison with Bayes factors

# Run bridge samplers
bridge_linear = bridge_sampler(fit_linear, silent=T)
bridge_quadratic = bridge_sampler(fit_quadratic, silent=T)

# Compute Bayes factors
bf(bridge_linear, bridge_quadratic)

