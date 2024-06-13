# # Making synthetic datasets: All use same variable ordering, just a question of how much of the training data used the 
# 
# # 1) Good synthesis: all of confidential data used
# # 2) Medium synthesis: 50% of minority confidential data used
# # 3) Bad synthesis: Only majority data used
# 
# library(tidyverse)
# library(ipumsr)
# library(srvyr)
# library(tidysynthesis)
# library(parsnip)
# library("data.table")
# library(recipes)
# library(devtools)
# load_all("../../syntheval")
# #library(syntheval)
# library(dials)
# library(tune)
# library(gtools)
# library(MASS)
# library(caret)
# library(gt)
# library(ggpubr)
# 
# set.seed(1)
# 
# # function that takes start data and confidential data, and synthesizes dataset
# synth <- function(start_data, conf_data, var_order){
#   # synthesize categorical first, then continuous
#   visit_sequence = visit_sequence(var_order, start_data, type = "manual")
#   roadmap = roadmap(conf_data, start_data, visit_sequence)
#   recipe = construct_recipes(roadmap = roadmap)
#   
#   tree_cl <- parsnip::decision_tree(cost_complexity = .0001) %>%
#     set_mode(mode = "classification") %>%
#     set_engine(engine = "rpart")
#   tree_reg <- parsnip::decision_tree(cost_complexity = .0001) %>%
#     set_mode(mode = "regression") %>%
#     set_engine(engine = "rpart")
#   
#   synth_algorithms = list()
#   for (i in 1:length(var_order)){
#     if(class(conf_data[,var_order[i]]) == "factor"){
#       synth_algorithms[[i]] = tree_cl
#     } else{
#       synth_algorithms[[i]] = tree_reg
#     }
#   }
#   
#   synth_spec = synth_spec(roadmap,
#                           synth_algorithms = synth_algorithms,
#                           recipe,
#                           predict_methods = sample_rpart)
#   
#   # noise
#   noise <- noise(roadmap = roadmap,
#                  add_noise = FALSE,
#                  exclusions = 0)
#   
#   # constraints
#   constraints <- constraints(roadmap = roadmap,
#                              constraints = NULL,
#                              max_z = 0)
#   
#   replicates <- replicates(replicates = 1,
#                            workers = 1,
#                            summary_function = NULL)
#   
#   # create a presynth object
#   presynth1 <- presynth(
#     roadmap = roadmap,
#     synth_spec = synth_spec,
#     noise = noise, 
#     constraints = constraints,
#     replicates = replicates
#   )
#   synthesized = synthesize(presynth1, progress = TRUE)
#   synthesized
# }
# 
# # load data
# data_mi <- read_csv("data/data_mi.csv")
# data_mi <- as.data.frame(unclass(data_mi),stringsAsFactors=TRUE)
# data_mi$pov <- as.factor(data_mi$pov)
# data_mi$BLACK <- as.factor(data_mi$BLACK)
# 
# # variable ordering: 
# dmy <- dummyVars(" ~ .", data = data_mi)
# trsf <- data.frame(predict(dmy, newdata = data_mi))
# cor <- data.frame(cor(trsf))
# rownames(cor)[order(abs(cor$BLACK.1), decreasing = TRUE)]
# 
# # based on this (select variables most highly correlated with BLACK first), and doing categorical variables first, order is:
# order = c("pov", "SEX", "EMPSTAT", "EDUC", "HCOVANY", "VETSTAT", "FTOTINC", "AGE")
# 
# ########## synthesize (Black, from all input data) ##########
# conf_data = data_mi
# start_data = data.frame("BLACK" = data_mi[,"BLACK"])
# synth_black_good <- synth(start_data, conf_data, order)
# 
# ########## synthesize (Black, from 50% of minority input data and all majority) ##########
# black_conf <- data_mi %>%
#   filter(BLACK == 1)
# # select random sample of these rows
# sample_black_conf <- sample_frac(black_conf, 0.5, replace = FALSE)
# not_black_conf = data_mi %>%
#   filter(BLACK == 0)
# conf_data <- rbind(sample_black_conf, not_black_conf)
# start_data = data.frame("BLACK" = as.factor(c(rep(0, nrow(not_black_conf)), 
#                                             rep(1, nrow(data_mi)-nrow(not_black_conf)))))
# synth_black_med <- synth(start_data, conf_data, order)
# 
# ########### synthesize (poverty, from only majority) ##########
# conf_data = data_mi %>%
#   filter(BLACK == 0)
# start_data = data.frame("BLACK" = as.factor(rep(0, nrow(data_mi))))
# synth_black_bad <- synth(start_data, conf_data, order)
# 
# # apply permutation to dataset (need to add Black people back in)
# synth_black_bad$synthetic_data$BLACK <- permute(data_mi$BLACK)
# 
# ########## synthesize (poverty, from 75% of minority input data and all majority) ##########
# black_conf <- data_mi %>%
#   filter(BLACK == 1)
# # select random sample of these rows
# sample_black_conf <- sample_frac(black_conf, 0.25, replace = FALSE)
# not_black_conf = data_mi %>%
#   filter(BLACK == 0)
# conf_data <- rbind(sample_black_conf, not_black_conf)
# start_data = data.frame("BLACK" = as.factor(c(rep(0, nrow(not_black_conf)), 
#                                             rep(1, nrow(data_mi)-nrow(not_black_conf)))))
# synth_black_75 <- synth(start_data, conf_data, order)
# 
# ########### save datasets ###########
# write_csv(synth_black_good$synthetic_data, "data/good_synth_black.csv")
# write_csv(synth_black_med$synthetic_data, "data/med_synth_black.csv")
# write_csv(synth_black_75$synthetic_data, "data/synth_75_black.csv")
# write_csv(synth_black_bad$synthetic_data, "data/bad_synth_black.csv")
