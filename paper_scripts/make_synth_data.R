# Making synthetic datasets: All use same variable ordering, just a question of how much of the training data used the 

# 1) Good synthesis: all of confidential data used
# 2) Medium synthesis: 50% of minority confidential data used
# 3) Bad synthesis: Only majority data used

library(tidyverse)
library(ipumsr)
library(srvyr)
library(tidysynthesis)
library(parsnip)
library("data.table")
library(recipes)
library(devtools)
load_all("../../syntheval")
#library(syntheval)
library(dials)
library(tune)
library(gtools)
library(MASS)
library(caret)
library(gt)
library(ggpubr)

set.seed(78483)

# function that takes start data and confidential data, and synthesizes dataset
synth <- function(start_data, conf_data, var_order){
  # synthesize categorical first, then continuous
  visit_sequence = visit_sequence(var_order, start_data, type = "manual")
  roadmap = roadmap(conf_data, start_data, visit_sequence)
  recipe = construct_recipes(roadmap = roadmap)
  
  tree_cl <- parsnip::decision_tree(cost_complexity = .0001) %>%
    set_mode(mode = "classification") %>%
    set_engine(engine = "rpart")
  tree_reg <- parsnip::decision_tree(cost_complexity = .0001) %>%
    set_mode(mode = "regression") %>%
    set_engine(engine = "rpart")
  
  synth_algorithms = list()
  for (i in 1:length(var_order)){
    if(class(conf_data[,var_order[i]]) == "factor"){
      synth_algorithms[[i]] = tree_cl
    } else{
      synth_algorithms[[i]] = tree_reg
    }
  }
  
  synth_spec = synth_spec(roadmap,
                          synth_algorithms = synth_algorithms,
                          recipe,
                          predict_methods = sample_rpart)
  
  # noise
  noise <- noise(roadmap = roadmap,
                 add_noise = FALSE,
                 exclusions = 0)
  
  # constraints
  constraints <- constraints(roadmap = roadmap,
                             constraints = NULL,
                             max_z = 0)
  
  replicates <- replicates(replicates = 1,
                           workers = 1,
                           summary_function = NULL)
  
  # create a presynth object
  presynth1 <- presynth(
    roadmap = roadmap,
    synth_spec = synth_spec,
    noise = noise, 
    constraints = constraints,
    replicates = replicates
  )
  synthesized = synthesize(presynth1, progress = TRUE)
  synthesized
}

# load data
data_mi <- read_csv("data/data_mi.csv")
data_mi <- as.data.frame(unclass(data_mi),stringsAsFactors=TRUE)
data_mi$pov <- as.factor(data_mi$pov)
data_mi$BLACK <- as.factor(data_mi$BLACK)

# variable ordering: 
dmy <- dummyVars(" ~ .", data = data_mi)
trsf <- data.frame(predict(dmy, newdata = data_mi))
cor <- data.frame(cor(trsf))
rownames(cor)[order(abs(cor$pov.1), decreasing = TRUE)]

# based on this (select variables most highly correlated with pov first), and doing categorical variables first, order is:
order = c("EMPSTAT", "EDUC", "BLACK", "SEX", "HCOVANY", "VETSTAT", "FTOTINC", "AGE")


permute_entries <- function(dt, k) {
  n <- nrow(dt)
  permuted_dt <- dt
  num_entries_to_permute <- round(n * k)
  
  for (col in names(dt)) {
    indices <- sample(1:n, num_entries_to_permute, replace = FALSE)
    indices_perm <- permute(indices)
    permuted_dt[indices_perm,col] <- permuted_dt[indices, col]
  }
  
  return(permuted_dt)
}

########## synthesize (poverty, from all input data) ##########
conf_data = data_mi
start_data = data.frame("pov" = data_mi[,"pov"])
synth_pov_good <- synth(start_data, conf_data, order)

######### synthesize (poverty, permute 50% of the entries in each column of poverty) ##########
pov_conf <- data_mi %>%
  filter(pov == 1)
# permute 50%
sample_pov_conf <- permute_entries(pov_conf, 0.25)
nopov_conf = data_mi %>%
  filter(pov == 0)
conf_data <- rbind(sample_pov_conf, nopov_conf)
start_data = data.frame("pov" = as.factor(c(rep(0, nrow(nopov_conf)), 
                                            rep(1, nrow(data_mi)-nrow(nopov_conf)))))
synth_pov_med <- synth(start_data, conf_data, order)

######### synthesize (poverty, permute 100% of the entries in each column of poverty) ##########
pov_conf <- data_mi %>%
  filter(pov == 1)
# permute 50%
sample_pov_conf <- permute_entries(pov_conf, .75)
nopov_conf = data_mi %>%
  filter(pov == 0)
conf_data <- rbind(sample_pov_conf, nopov_conf)
start_data = data.frame("pov" = as.factor(c(rep(0, nrow(nopov_conf)), 
                                            rep(1, nrow(data_mi)-nrow(nopov_conf)))))
synth_pov_bad <- synth(start_data, conf_data, order)

######### synthesize (poverty, permute 50% of the entries in each column of poverty) ##########
pov_conf <- data_mi %>%
  filter(pov == 1)
# permute 50%
nopov_conf = data_mi %>%
  filter(pov == 0)
sample_nopov_conf <- permute_entries(nopov_conf, 0.25)
conf_data <- rbind(sample_nopov_conf, pov_conf)
start_data = data.frame("pov" = as.factor(c(rep(0, nrow(nopov_conf)), 
                                            rep(1, nrow(data_mi)-nrow(nopov_conf)))))
synth_nopov_med <- synth(start_data, conf_data, order)

######### synthesize (poverty, permute 100% of the entries in each column of no poverty) ##########
pov_conf <- data_mi %>%
  filter(pov == 1)
# permute 50%
nopov_conf = data_mi %>%
  filter(pov == 0)
sample_nopov_conf <- permute_entries(nopov_conf, .75)
conf_data <- rbind(sample_nopov_conf, pov_conf)
start_data = data.frame("pov" = as.factor(c(rep(0, nrow(nopov_conf)), 
                                            rep(1, nrow(data_mi)-nrow(nopov_conf)))))
synth_nopov_bad <- synth(start_data, conf_data, order)

write_csv(synth_pov_good$synthetic_data, "data/good_synth_pov2.csv")
write_csv(synth_pov_med$synthetic_data, "data/med_synth_pov2.csv")
write_csv(synth_pov_bad$synthetic_data, "data/bad_synth_pov2.csv")
write_csv(synth_nopov_med$synthetic_data, "data/med_synth_nopov2.csv")
write_csv(synth_nopov_bad$synthetic_data, "data/bad_synth_nopov2.csv")





















########## synthesize (poverty, from 50% of minority input data and all majority) ##########
pov_conf <- data_mi %>%
  filter(pov == 1)
# select random sample of these rows
sample_pov_conf <- sample_frac(pov_conf, 0.5, replace = FALSE)
nopov_conf = data_mi %>%
  filter(pov == 0)
conf_data <- rbind(sample_pov_conf, nopov_conf)
start_data = data.frame("pov" = as.factor(c(rep(0, nrow(nopov_conf)), 
                                            rep(1, nrow(data_mi)-nrow(nopov_conf)))))
synth_pov_med <- synth(start_data, conf_data, order)

########### synthesize (poverty, from only majority) ##########
conf_data = data_mi %>%
  filter(pov == 0)
start_data = data.frame("pov" = as.factor(rep(0, nrow(data_mi))))
synth_pov_bad <- synth(start_data, conf_data, order)

# apply permutation to dataset (need to add people in poverty back in)
synth_pov_bad$synthetic_data$pov <- permute(data_mi$pov)

########## synthesize (poverty, from 75% of minority input data and all majority) ##########
pov_conf <- data_mi %>%
  filter(pov == 1)
# select random sample of these rows
sample_pov_conf <- sample_frac(pov_conf, 0.25, replace = FALSE)
nopov_conf = data_mi %>%
  filter(pov == 0)
conf_data <- rbind(sample_pov_conf, nopov_conf)
start_data = data.frame("pov" = as.factor(c(rep(0, nrow(nopov_conf)), 
                                            rep(1, nrow(data_mi)-nrow(nopov_conf)))))
synth_pov_75 <- synth(start_data, conf_data, order)

########## synthesize (poverty, from 75% of minority input data and all majority) ##########
pov_conf <- data_mi %>%
  filter(pov == 1)
# select random sample of these rows
sample_pov_conf <- sample_frac(pov_conf, 0.75, replace = FALSE)
nopov_conf = data_mi %>%
  filter(pov == 0)
conf_data <- rbind(sample_pov_conf, nopov_conf)
start_data = data.frame("pov" = as.factor(c(rep(0, nrow(nopov_conf)), 
                                            rep(1, nrow(data_mi)-nrow(nopov_conf)))))
synth_pov_25 <- synth(start_data, conf_data, order)

########### save datasets ###########
write_csv(synth_pov_good$synthetic_data, "data/good_synth_pov.csv")
write_csv(synth_pov_med$synthetic_data, "data/med_synth_pov.csv")
write_csv(synth_pov_75$synthetic_data, "data/synth_pov_75.csv")
write_csv(synth_pov_25$synthetic_data, "data/synth_pov_25.csv")
write_csv(synth_pov_bad$synthetic_data, "data/bad_synth_pov.csv")
