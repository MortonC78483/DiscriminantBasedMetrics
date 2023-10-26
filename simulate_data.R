## ---------------------------
##
## Script name: simulate_data
##
## Purpose of script: Create simulated datasets for confidential and synthetic data
##
## Author: Claire Morton
##
## Date Created: 2032-10-24
##
## Email: mortonc@stanford.edu
##
## ---------------------------
##
## Notes:
##   I will simulate datasets, first with 10 variables and 10,000 data points, 
##   that contain only continuous, then only (binary) categorical, data points. 
##   I will create several versions of simulated confidential/synthetic data to 
##   represent different strengths and weaknesses that could occur in a synthesized dataset.
##   
##  1) Confidential and synthetic dataset are drawn from (for continuous) multivariate 
##    normal with known variance-covariance matrix, and from (for categorical) known 
##    distribution. Control pairing to test which hyperparameters work best in continuous 
##    vs. categorical “ideal” setting.
##  2) Confidential data are drawn from (for continuous) multivariate normal with known 
##    variance-covariance matrix, and from (for categorical) known distribution. 
##    Synthetic data are drawn from distribution with larger variances and slightly (~20%) 
##    off means.
##
## ---------------------------
library(MASS)
library(GenOrd)

#==========================================================================
set.seed(101)
n=10000 # number of samples
corr=.7 # correlation
p = 10 # number of variables

#==========================================================================
# setting 1 (control)
# continuous
vcov = matrix(corr, p, p) + diag(1-corr, p, p) # vcov matrix, p on all off-diagonal variables
means = rep(0,p) # mean vector
conf_cont1 = mvrnorm(means, vcov, n = n)
synth_cont1 = mvrnorm(means, vcov, n = n)

# categorical
marginal = rep(list(c(.5)), p)
corrcheck(marginal) # check ok
vcov = vcov
conf_cat1 = (2*ordsample(n, marginal, vcov))-3 # -1 and 1
synth_cat1 = (2*ordsample(n, marginal, vcov))-3 # -1 and 1

#==========================================================================
# setting 2 (everything slightly off)
# continuous
conf_cont2 = conf_cont1
vcov_synth = vcov # increase variance
means_synth = jitter(means, amount = .4) # slightly offset means
synth_cont2 = mvrnorm(means_synth, vcov_synth, n = n)

# categorical
conf_cat2 = conf_cat1
vcov_cat = vcov/2+diag(1/2, 10, 10)
marginal = sapply((1-means_synth)/2, list) # set mean to be offset means above
corrcheck(marginal) # check ok
synth_cat2 = (2*ordsample(n, marginal, vcov_cat))-3 # -1 and 1

#==========================================================================
# logistic regression (not tidymodels)
# setting 1, logistic regression
combined_1 = data.frame(rbind(cbind(conf_cont1, "confidential" = 1),
                              cbind(synth_cont1, "confidential" = 0)))
combined_1$confidential = as.factor(combined_1$confidential)
model1 <- glm(confidential ~.,family=binomial(link='logit'),data=combined_1)
prop_scores1 = model1$fitted.values


# setting 2, logistic regression
combined_2 = data.frame(rbind(cbind(conf_cont2, "confidential" = 1),
                              cbind(synth_cont2, "confidential" = 0)))
combined_2$confidential = as.factor(combined_2$confidential)
model2 <- glm(confidential ~.,family=binomial(link='logit'),data=combined_2)
prop_scores2 = model2$fitted.values

# logistic regression (tidymodels)

# CART (not tidymodels)

# CART(tidymodels)

#==========================================================================
# ggplot(data.frame(prop_scores1), aes(x = prop_scores1))+
#   geom_density()+
#   geom_density(data = data.frame(prop_scores2), aes(x = prop_scores2), color = "red")
# 
# hist(combined_2$V10)
# hist(prop_scores2)

# example
# data = data.frame(rbind(cbind(rep(1, 1000), "confidential"=1),
#                         cbind(rep(0, 1000), "confidential"=0)))
# data[1,] = c(0, 1)
# data[1001,] = c(1,0)
# model = glm(confidential ~.,family=binomial(link='logit'),data=data)
# prop_scores = model$fitted.values
