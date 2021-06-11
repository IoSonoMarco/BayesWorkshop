library(plyr)
library(dplyr)
library(brms)
library(ggplot2)
library(tidybayes)
library(modelr) 
library(ggdist)
library(bridgesampling)
library(loo)

load('fits/bayesian_hierarchical_changepoint_loos.Rda')
load('fits/bayesian_hierarchical_changepoint_models.Rda')

### load and visualize data ###
data <- read.csv('datasets/CP.csv')

ggplot(data, aes(x=Numerosity, y=RT, group=Subj)) +
  geom_point() 

# b1: pre-inflection slope
# b2: post-inflection slope
bform <- brms::bf( # explicit brms dependencies (bf is also a function of bridgesampling)
  RT ~ b0 + b1 * (Numerosity - 5) * step(5 - Numerosity) + b2 * (Numerosity - 5) * step(Numerosity - 5),
  b0 + b1 + b2 ~ 1 + (1|Subj),
  nl = TRUE
)

bprior <- set_prior("normal(0, 100)", lb = 0, nlpar = 'b0') +
          prior(normal(0, 100), nlpar = "b1") +
          prior(normal(0, 100), nlpar = "b2") 
          
m <- brm(formula = bform, data = data, prior = bprior, iter = 4000, warmup = 500, chains = 3, cores = 3)

data$pred_m0 <- fitted(m)

ggplot(data, aes(x=Numerosity, y=pred_m0[,'Estimate'], group=Subj)) +
  geom_line()

data %>%
  add_predicted_draws(m, re_formula = NA) %>%
  ggplot(aes(x = Numerosity, y = RT)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .90, .8, .5), color = "blue") +
  geom_point(data = data, size = 2, alpha=0.5) +
  scale_fill_brewer() +
  theme_bw()


### Group difference analysis ###

ggplot(data, aes(x=Numerosity, y=RT, group=Subj, color=Group)) +
  facet_grid(~Group, labeller = label_both) +
  geom_point() 

### Model with changepoint at Numerosity=3 ###

bform <- brms::bf( 
  RT ~ b0 + b1 * (Numerosity - 3) * step(3 - Numerosity) + b2 * (Numerosity - 3) * step(Numerosity - 3),
  b0 + b1 + b2 ~ 1 + Group + (1|Subj),
  nl = TRUE
)

bprior <- set_prior("normal(0, 100)", lb = 0, nlpar = 'b0') +
  prior(normal(0, 100), nlpar = "b1") +
  prior(normal(0, 100), nlpar = "b2") 

m_cp3 <- brm(formula = bform, data = data, prior = bprior, iter = 4000, warmup = 500, chains = 3, cores = 3, 
             save_pars=save_pars(all = TRUE), sample_prior = TRUE)

### Model with changepoint at Numerosity=4 ###

bform <- brms::bf( 
  RT ~ b0 + b1 * (Numerosity - 4) * step(4 - Numerosity) + b2 * (Numerosity - 4) * step(Numerosity - 4),
  b0 + b1 + b2 ~ 1 + Group + (1|Subj),
  nl = TRUE
)

bprior <- set_prior("normal(0, 100)", lb = 0, nlpar = 'b0') +
  prior(normal(0, 100), nlpar = "b1") +
  prior(normal(0, 100), nlpar = "b2") 

m_cp4 <- brm(formula = bform, data = data, prior = bprior, iter = 4000, warmup = 500, chains = 3, cores = 3, 
             save_pars=save_pars(all = TRUE), sample_prior = TRUE)

### Model with changepoint at Numerosity=6 ###

bform <- brms::bf( 
  RT ~ b0 + b1 * (Numerosity - 6) * step(6 - Numerosity) + b2 * (Numerosity - 6) * step(Numerosity - 6),
  b0 + b1 + b2 ~ 1 + Group + (1|Subj),
  nl = TRUE
)

bprior <- set_prior("normal(0, 100)", lb = 0, nlpar = 'b0') +
  prior(normal(0, 100), nlpar = "b1") +
  prior(normal(0, 100), nlpar = "b2") 

m_cp6 <- brm(formula = bform, data = data, prior = bprior, iter = 4000, warmup = 500, chains = 3, cores = 3, 
             save_pars=save_pars(all = TRUE), sample_prior = TRUE)

# leave-one-out

loo_cp3 <- loo(m_cp3) # reloo = TRUE
loo_cp4 <- loo(m_cp4)
loo_cp6 <- loo(m_cp6)

# Inspect comparison
loo_compare(loo_cp3, loo_cp4, loo_cp6)

#####################################
### analyze m_cp4, the best model ###
#####################################

### posterior predictive check ###

pp_check(m_cp4)
pp_check(m_cp4, type='stat', stat='mean')

post_pred <- add_predicted_draws(data, m_cp4) %>%
  ddply(.(Group, Subj, Numerosity), summarize, post_mean_RT=mean(.prediction))

ggplot(post_pred, aes(x=Numerosity, y=post_mean_RT, group=Subj, color=Group)) +
  facet_grid(~ Group, labeller = label_both) +
  geom_point(alpha=0.8) +
  geom_point(data=data, mapping=aes(x=Numerosity, y=RT, group=Subj), color='black', alpha=0.3) 

### contrast analysis ###

# set the hypotheses #

contr <- c('pre-slope G vs pre-slope A' = 'b1_GroupG = 0',
           'post-slope G vs post-slope A' = 'b2_GroupG = 0',
           'pre-slope A vs post-slope A' = 'b2_Intercept - b1_Intercept = 0',
           'pre-slope G vs post-slope G' = '(b2_Intercept + b2_GroupG) - (b1_Intercept + b1_GroupG) = 0',
           'diff A - B' = '((b2_Intercept + b2_GroupG) - (b1_Intercept + b1_GroupG)) = (b2_Intercept - b1_Intercept)')
h <- hypothesis(m_cp4, contr)
plot(h)[[1]] + geom_vline(xintercept = 0, col='red', linetype='dashed')

#####################################
### analyze m_cp4, the best model ###
#####################################

bform <- brms::bf( 
  RT ~ b0 + b1 * (Numerosity - 4) * step(4 - Numerosity) + b2 * (Numerosity - 4) * step(Numerosity - 4),
  b0 + b1 + b2 ~ 1 + Group + (1|gr(Subj, by=Group)),
  nl = TRUE
)

bprior <- set_prior("normal(0, 100)", lb = 0, nlpar = 'b0') +
  prior(normal(0, 100), nlpar = "b1") +
  prior(normal(0, 100), nlpar = "b2") 

m_cp4_cond_cov <- brm(formula = bform, data = data, prior = bprior, iter = 4000, warmup = 500, chains = 3, cores = 3, 
             save_pars=save_pars(all = TRUE), sample_prior = TRUE)

### model selection (parsimony) ###

bridge_0 <- bridge_sampler(m_cp4)
bridge_1 <- bridge_sampler(m_cp4_cond_cov)
bf(bridge_0, bridge_1)

