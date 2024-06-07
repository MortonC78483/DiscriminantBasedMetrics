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
order = c("pov", "SEX", "EMPSTAT", "EDUC", "HCOVANY", "VETSTAT", "FTOTINC", "AGE")


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
start_data = data.frame("BLACK" = data_mi[,"BLACK"])
synth_BLACK_good <- synth(start_data, conf_data, order)

######### synthesize (poverty, permute 50% of the entries in each column of poverty) ##########
BLACK_conf <- data_mi %>%
  filter(BLACK == 1)
# permute 50%
sample_BLACK_conf <- permute_entries(BLACK_conf, 0.25)
noBLACK_conf = data_mi %>%
  filter(BLACK == 0)
conf_data <- rbind(sample_BLACK_conf, noBLACK_conf)
start_data = data.frame("BLACK" = as.factor(c(rep(0, nrow(noBLACK_conf)), 
                                            rep(1, nrow(data_mi)-nrow(noBLACK_conf)))))
synth_BLACK_med <- synth(start_data, conf_data, order)

######### synthesize (BLACKerty, permute 100% of the entries in each column of BLACKerty) ##########
BLACK_conf <- data_mi %>%
  filter(BLACK == 1)
# permute 50%
sample_BLACK_conf <- permute_entries(BLACK_conf, .75)
noBLACK_conf = data_mi %>%
  filter(BLACK == 0)
conf_data <- rbind(sample_BLACK_conf, noBLACK_conf)
start_data = data.frame("BLACK" = as.factor(c(rep(0, nrow(noBLACK_conf)), 
                                            rep(1, nrow(data_mi)-nrow(noBLACK_conf)))))
synth_BLACK_bad <- synth(start_data, conf_data, order)

######### synthesize (BLACKerty, permute 50% of the entries in each column of BLACKerty) ##########
BLACK_conf <- data_mi %>%
  filter(BLACK == 1)
# permute 50%
noBLACK_conf = data_mi %>%
  filter(BLACK == 0)
sample_noBLACK_conf <- permute_entries(noBLACK_conf, 0.25)
conf_data <- rbind(sample_noBLACK_conf, BLACK_conf)
start_data = data.frame("BLACK" = as.factor(c(rep(0, nrow(noBLACK_conf)), 
                                            rep(1, nrow(data_mi)-nrow(noBLACK_conf)))))
synth_noBLACK_med <- synth(start_data, conf_data, order)

######### synthesize (BLACKerty, permute 100% of the entries in each column of no BLACKerty) ##########
BLACK_conf <- data_mi %>%
  filter(BLACK == 1)
# permute 50%
noBLACK_conf = data_mi %>%
  filter(BLACK == 0)
sample_noBLACK_conf <- permute_entries(noBLACK_conf, .75)
conf_data <- rbind(sample_noBLACK_conf, BLACK_conf)
start_data = data.frame("BLACK" = as.factor(c(rep(0, nrow(noBLACK_conf)), 
                                            rep(1, nrow(data_mi)-nrow(noBLACK_conf)))))
synth_noBLACK_bad <- synth(start_data, conf_data, order)

write_csv(synth_BLACK_good$synthetic_data, "data/good_synth_BLACK2.csv")
write_csv(synth_BLACK_med$synthetic_data, "data/med_synth_BLACK2.csv")
write_csv(synth_BLACK_bad$synthetic_data, "data/bad_synth_BLACK2.csv")
write_csv(synth_noBLACK_med$synthetic_data, "data/med_synth_noBLACK2.csv")
write_csv(synth_noBLACK_bad$synthetic_data, "data/bad_synth_noBLACK2.csv")
