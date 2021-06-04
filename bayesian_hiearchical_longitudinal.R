library(plyr)
library(dplyr)
library(brms)
library(ggplot2)
library(bayesplot)
library(bridgesampling)
library(tidybayes)
library(bayestestR)

#load('fits/bayesian_hierarchical_longitudinal_m1.Rda')
#load('fits/bayesian_hierarchical_longitudinal_m2.Rda')
#load('fits/bayesian_hierarchical_longitudinal_m3.Rda')
#load('fits/bayesian_hierarchical_longitudinal_m4.Rda')

### load and visualize data ###
data <- read.csv('datasets/IGT_Block.csv') %>%
  mutate(Subj=factor(Subj), Good_prop = Good/Total, Block=Block-1) %>%
  filter(Condition != '1') %>%
  mutate(Condition=factor(Condition))
  

data_stats <- ddply(data, .(Condition, Block), summarize, mean_prop=mean(Good_prop), se=sd(Good_prop)/sqrt(length(Good_prop)))

ggplot(data_stats, aes(x=Block, y=mean_prop, group=Condition, color=Condition)) +
  geom_pointrange(aes(ymin=mean_prop-se, ymax=mean_prop+se)) +
  geom_line()

### specify models ###

###########################
### parallel sub-models ###
###########################

m1 <- brm(Good|trials(Total) ~ Block + Condition + (1|Subj), 
          family=binomial(link='logit'), data=data, 
          iter=3500, warmup=500, chains=3, cores=3,
          save_pars=save_pars(all = TRUE))

print(m1)

draws <- as.matrix(m1)
mcmc_areas(draws, pars=c('b_Intercept', 'b_Block', 'b_Condition3'), prob=0.95)

data$pred_m1 <- fitted(m1, scale = 'response', re_formula = NA) 

### plot predicted data ###

ggplot(data, aes(x=Block, y=pred_m1[,'Estimate'], group=Condition, color=Condition)) +
  geom_pointrange(aes(ymin=pred_m1[,'Q2.5'], ymax=pred_m1[,'Q97.5']), position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2))

###############################
### non-parallel sub-models ###
###############################

m2 <- brm(Good|trials(Total) ~ Block*Condition + (1 + Block|Subj), 
          family=binomial(link='logit'), data=data, 
          iter=3500, warmup=500, chains=3, cores=3,
          save_pars=save_pars(all = TRUE))

### compare m1 and m2 via Bayes Factor ###
bridge <- bridge_sampler(m1)

bridge2 <- bridge_sampler(m2)

bf(bridge2, bridge)

###################################
### conditional quadratic model ###
###################################

### add squared predictor ###
data$Block_sq <- data$Block^2

m3 <- brm(Good|trials(Total) ~ Block + Block_sq*Condition + (1 + Block + Block_sq|Subj), 
          family=binomial(link='logit'), data=data, 
          iter=3500, warmup=500, chains=3, cores=3,
          save_pars=save_pars(all = TRUE))

bridge3 <- bridge_sampler(m3)

bf(bridge3, bridge2)

################################
### complete quadratic model ###
################################

m4 <- brm(Good|trials(Total) ~ (Block + Block_sq)*Condition + (1 + Block + Block_sq|Subj), 
          family=binomial(link='logit'), data=data, 
          iter=3500, warmup=500, chains=3, cores=3,
          save_pars=save_pars(all = TRUE))

bridge4 <- bridge_sampler(m4)

bf(bridge4, bridge3)

barplot(post_prob(bridge, bridge2, bridge3, bridge4))

### after model selection ###

print(m4)

conditional_effects(m4)

data$pred_m4 <- fitted(m4, scale = 'response') #, re_formula = NA)

ggplot(data, aes(x=Block, y=pred_m4[,'Estimate'], group=Condition, color=Condition)) +
  geom_pointrange(aes(ymin=pred_m4[,'Q2.5'], ymax=pred_m4[,'Q97.5']), position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2))

### predict data for contrast analysis ###

# with RE #
pp <- add_predicted_draws(data, m4)

pp_stats <- ddply(pp, .(Condition, Subj, Block), summarize, m=mean(.prediction))

ggplot(pp_stats, aes(x=Block, y=m, group=Subj, col=Condition)) +
  facet_grid(~ Condition) +
  geom_point()

# without RE (only FE) #
pp <- add_predicted_draws(data, m4, re_formula = NA)

pp_stats <- ddply(pp, .(Condition, Subj, Block), summarize, m=mean(.prediction))

ggplot(pp_stats, aes(x=Block, y=m, group=Subj, col=Condition)) +
  facet_grid(~ Condition) +
  geom_point()

## contrast: Condition 2 VS Condition 3 | Block 2 ##

pp <- add_predicted_draws(data, m4, re_formula = NA)

# plot predictions #
pp_stats <- ddply(pp, .(Condition, Block), summarize, mean=mean(.prediction), ci_low=ci(.prediction)$CI_low, ci_high=ci(.prediction)$CI_high)

ggplot(pp_stats, aes(x=Block, y=mean, group=Condition, color=Condition)) +
  geom_point(position=position_dodge(0.2), size=3) +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=0, position=position_dodge(0.2), size=2, alpha=0.3) +
  geom_line(position=position_dodge(0.2))

# do contrasts on predictions #

block_2_pred <- pp[pp$Block==2,]
cond_2 <- block_2_pred$.prediction[block_2_pred$Condition=='2']
cond_3 <- block_2_pred$.prediction[block_2_pred$Condition=='3']

hist(cond_2)
hist(cond_3)

hist(cond_2 - cond_3)

cond_2_pred <- pp[pp$Condition=='2',]
block_0 <- cond_2_pred$.prediction[cond_2_pred$Block==0]
block_3 <- cond_2_pred$.prediction[cond_2_pred$Block==3]

hist(block_0)
hist(block_3)

hist(block_0 - block_3)

ci(block_0 - block_3, method='ETI')

# with posterior_predict #

b0 <- posterior_predict(m4, data[data$Condition=='2' & data$Block==0,], re_formula = NA)
b3 <- posterior_predict(m4, data[data$Condition=='2' & data$Block==3,], re_formula = NA)
hist(c(b0)-c(b3))

