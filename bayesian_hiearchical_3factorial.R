library(dplyr)
library(brms)
library(ggplot2)
library(emmeans)

### load and visualize data ###
data <- read.csv('datasets/3FactRT.csv') 

# reduce sample size for didactic purpose #
data <- data[data$Subj_ID %in% seq(1,5),] %>%
  mutate(Polarity=factor(Polarity), Truth_Value=factor(Truth_Value), Numerosity=factor(Numerosity), Subj_ID=factor(Subj_ID))

# plot data group-level #
ggplot(data, aes(x=RT, col=Numerosity, fill=Numerosity)) +
  facet_grid(Polarity ~ Truth_Value, labeller = label_both) +
  geom_density(alpha=0.2) +
  xlim(0,3000)

# plot data individual-level #
ggplot(data, aes(x=Subj_ID, y=RT, fill=Numerosity)) +
  facet_grid(Polarity ~ Truth_Value, labeller = label_both) +
  geom_boxplot(alpha=0.5) +
  ylim(0,3000)

### transform dependent variable ###
data$RT_log <- log(data$RT)

ggplot(data, aes(x=RT_log, col=Numerosity, fill=Numerosity)) +
  facet_grid(Polarity ~ Truth_Value, labeller = label_both) +
  geom_density(alpha=0.2) 

#
get_prior(RT_log ~ Truth_Value*Numerosity + Polarity + (1|Subj_ID), data)
# 

### specify priors ###
priors <- prior(normal(0,100), class = "b") +
          prior(normal(0,100), class='Intercept')

### fit the model ###
m <- brm(RT_log ~ Truth_Value*Numerosity + Polarity + (1|Subj_ID), prior = priors,
         data = data, iter = 4000, warmup = 500, chains = 3, core = 3, 
         sample_prior = TRUE, control=list(adapt_delta=0.95))

#
prior_summary(m)
#

print(m)

plot(m)

### stats on posteriors ###
emm <- emmeans(m, specs = pairwise ~ Truth_Value*Numerosity)
print(emm$contrasts)
plot(emm$contrasts) + geom_vline(xintercept=0, col='red', linetype='dashed')

### contrast analysis on parameters ###

# extract posterior samples #
draws <- prepare_predictions(m)
str(draws)
# extract regression coefficients #
fe_draws <- draws$dpars$mu$fe$b

# manually extract posterior samples for contrast #
draws_True_Mag <- fe_draws[, 'b_Intercept'] + fe_draws[, 'b_Truth_ValueTrue']
draws_True_Min <- fe_draws[, 'b_Intercept'] + fe_draws[, 'b_Truth_ValueTrue'] + fe_draws[,'b_NumerosityMIN'] + fe_draws[,'b_Truth_ValueTrue:NumerosityMIN']
plot(density(draws_True_Mag), main='Posteriors', xlab='coeff')
lines(density(draws_True_Min), col='red')
# compute contrast #
contr_True_Mag_vs_Min <- draws_True_Mag - draTrue_Mws_in
plot(density(contr_True_Mag_vs_Min))
median(contr_True_Mag_vs_Min)
emm$contrasts # compare

### use hypothesis ###
contr <- c('True MAG - True MIN' = '- (NumerosityMIN + Truth_ValueTrue:NumerosityMIN) = 0')
H <- hypothesis(m, c(contr))
plot(H)[[1]] + xlim(-1,1)



h1 <- hypothesis(m, 'Truth_ValueTrue = 0')
plot(h1)[[1]] + xlim(-1,1)

h2 <- hypothesis(m, 'PolarityNEG = 0')
plot(h2)[[1]] + xlim(-1,1)

hh <- hypothesis(m, c('Truth_ValueTrue = 0',
                      'PolarityNEG = 0'))
