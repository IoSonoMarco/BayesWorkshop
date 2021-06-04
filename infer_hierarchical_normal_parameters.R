library(tidyr)
library(rstan)
library(bayesplot)
library(ggplot2)


### Simulate hierarchical data
set.seed(100)
group_n <- c(50, 3, 45, 40, 30, 25, 70, 60, 65, 75, 18) # observations for each group
N <- sum(group_n) # total number of data points
K <- length(group_n) # number of groups
mu <- 0
tau <- 3
sigma <- 2
group_betas <- rnorm(K, mu, tau) # sample betas (group-specific means) from the overarching normal distribution

# generate data fro each group
x <- c()
for(k in 1:K){ # for each group
  group_data <- rnorm(group_n[k], group_betas[k], sigma) # sample group data from a normal distribution
  x <- c(x, group_data) # stack the group data horizontally
}

# generate group indexes
ind <- c()
for(k in 1:K){ # for each group
  group_index <- rep(k, group_n[k]) # generate the index vector for the specific group
  ind <- c(ind, group_index) # stack indexes horizontally
}

d <- data.frame(group = factor(ind), score = x)
boxplot(score ~ group, data=d)

##################
### NO POOLING ###
##################

### Compile models
model_nopooling = stan_model('models/model_hierarchical_normal_nopooling.stan')

### Prepare a list of the data
data = list(
  x = x,
  I = ind,
  N = N,
  K = K
)

### Sample from models
fit_nopooling = sampling(model_nopooling, data=data, chains=2, iter=5000, warmup=500)

### Inspect summaries and R_hat diagnostics
print(fit_nopooling)

### visualization
draws <- as.matrix(fit_nopooling)
mcmc_areas(draws, pars=c(paste('beta[', seq(1,K), ']', sep='')), prob=0.95)

draws_long_nopooling <- draws %>%
  as.data.frame() %>%
  gather(parameter, samples, 'beta[1]':'beta[11]', factor_key=TRUE) 

ggplot(draws_long_nopooling[1:9000*4,], aes(x=samples, fill=parameter)) +
  geom_density(alpha=0.3) 

#######################
### PARTIAL POOLING ###
#######################

### Compile models
model_partialpooling = stan_model('models/model_hierarchical_normal_partialpooling.stan')

### Prepare a list of the data
data = list(
  x = x,
  I = ind,
  N = N,
  K = K
)

### Sample from models
fit_partialpooling = sampling(model_partialpooling, data=data, chains=2, iter=5000, warmup=500)

### Inspect summaries and R_hat diagnostics
print(fit_partialpooling)

### visualization
draws <- as.matrix(fit_partialpooling)
mcmc_areas(draws, pars=c(paste('beta[', seq(1,K), ']', sep='')), prob=0.95)

mcmc_areas(draws, pars=c('mu', 'tau'), prob=0.95)

mcmc_pairs(fit_partialpooling, pars=c(paste('beta[', seq(1,K), ']', sep='')), diag_fun="dens")

draws_long_partialpooling <- draws %>%
  as.data.frame() %>%
  gather(parameter, samples, 'beta[1]':'beta[11]', factor_key=TRUE)


library(ggdist)
ggplot(draws_long_nopooling[1:9000*4,], aes(x=samples, y=parameter)) +
  stat_halfeye(aes(fill='blue'), alpha=0.4) + 
  stat_halfeye(data=draws_long_partialpooling[1:9000*4,], aes(fill='green'), alpha=0.4) +
  scale_fill_manual(name='Legend', labels = c("No pooling", "Partial pooling"), values = c("blue", "green")) +
  labs(x="Parameter value", y='Parameter name') +
  geom_vline(xintercept = 0.06)


