# /*************************/
# Programmer: Claire Morton
# Date created: 06/12/2024
# Date of last revision: 06/12/2024
# original data: 
# 
# Description: Contains functions written for the project pulled by the other files
# 
# /*************************/
  
  ################################ GENERAL ################################



################################ MODEL COMPARISON ################################
# function creates a postsynth object from a dataset
# data: the input dataset that is returned as a postsynth object
generate_postsynth <- function(data){
  list(
    synthetic_data = data,
    jth_synthesis_time = data.frame(
      variable = factor(colnames(data.frame(data)))
    )
  )  %>%
    structure(class = "postsynth")
}

# function creates a confidential/synthetic dataset where the confidential and synthetic data are drawn from the same underlying
# distribution.
# the confidential dataset is a multivariate normal distribution
# p: the number of variables in the synthetic/confidental dataset
# n: the number of rows in the synthetic/confidental dataset
# corr: the value of the variance/covariance matrix on off-diagonals (1 on diagonals)
# return: list containing a confidential and a synthetic (as post synth object) dataset
generate_controls <- function(p, corr, n){
  # Same distribution
  vcov = matrix(corr, p, p) + diag(1-corr, p, p) # vcov matrix, p on all off-diagonal variables
  means = rep(0,p) # mean vector
  conf1 = as_tibble(mvrnorm(means, vcov, n = n), n = p)
  synth1 = generate_postsynth(as_tibble(mvrnorm(means, vcov, n = n), n=p))
  list(conf1, synth1)
}

# function creates a confidential/synthetic dataset where the synthetic data have slightly different means and 
# larger variances/covariances than the confidential dataset.
# p: the number of variables in the synthetic/confidental dataset
# n: the number of rows in the synthetic/confidental dataset
# corr: the value of the variance/covariance matrix on off-diagonals (1 on diagonals) for the confidential dataset.
#       the variance/covariance matrix of the synthetic data has 4 on the diagonals and 2*corr on the off-diagonals
# prop: the proportion of the synthetic dataset that should not be drawn from the same distribution as the confidential dataset
# return: list containing a confidential and a synthetic (as post synth object) dataset
generate_alloff <- function(p, corr, n, prop = 1){
  # larger variances, slightly off means
  vcov = matrix(corr, p, p) + diag(1-corr, p, p) # vcov matrix, p on all off-diagonal variables
  means = rep(0,p) # mean vector
  conf2 = as_tibble(mvrnorm(means, vcov, n = n), n=p)
  means_new = means
  means_new = sample(c(-.2,.2), 10, replace = T) # mean vector
  vcov_new = vcov
  vcov_new = matrix(2*corr, p, p) + diag(4-2*corr, p, p)
  if (prop != 1){
    synth2 = generate_postsynth(rbind(sample_n(as_tibble(mvrnorm(means_new, vcov_new, n = n), n=p), size=n*prop, replace = FALSE),
                                      as_tibble(mvrnorm(means, vcov, n = n*(1-prop)), n=p)))
  } else{
    synth2 = generate_postsynth(rbind(sample_n(as_tibble(mvrnorm(means_new, vcov_new, n = n), n=p), size=n*prop, replace = FALSE)))
  }
  
  list(conf2, synth2)
}

# function creates a confidential/synthetic dataset where the first variable in the synthetic dataset has a different
# mean than the corresponding variable in the confidential dataset (2 instead of 0)
# p: the number of variables in the synthetic/confidental dataset
# n: the number of rows in the synthetic/confidental dataset
# corr: the value of the variance/covariance matrix on off-diagonals (1 on diagonals) for the confidential and synthetic datasets
# prop: the proportion of the synthetic dataset that should not be drawn from the same distribution as the confidential dataset
# return: list containing a confidential and a synthetic (as post synth object) dataset
generate_mean <- function(p, corr, n, prop = 1){
  # one variable centered at wrong value
  vcov = matrix(corr, p, p) + diag(1-corr, p, p) # vcov matrix, p on all off-diagonal variables
  means = rep(0,p) # mean vector
  conf3 = as_tibble(mvrnorm(means, vcov, n = n), n=p)
  means_new = means
  means_new[1] = 2
  if (prop != 1){
    synth3 = generate_postsynth(rbind(sample_n(as_tibble(mvrnorm(means_new, vcov, n = n), n=p), size=n*prop, replace = FALSE),
                                      as_tibble(mvrnorm(means, vcov, n = n*(1-prop)), n=p)))
  } else{
    synth3 = generate_postsynth(rbind(sample_n(as_tibble(mvrnorm(means_new, vcov, n = n), n=p), size=n*prop, replace = FALSE)))
  }
  
  list(conf3, synth3)
}

# function creates a confidential/synthetic dataset where the first variable in the synthetic dataset has a negative, rather than
# positive, covariance with the other variables in the dataset
# p: the number of variables in the synthetic/confidental dataset
# n: the number of rows in the ynthetic/confidental dataset
# corr: the value of the variance/covariance matrix on off-diagonals (1 on diagonals) for the confidential and synthetic datasets
#       this value is negated in the simulated synthetic dataset's covariance matrix for the first row and column
# prop: the proportion of the synthetic dataset that should not be drawn from the same distribution as the confidential dataset
# return: list containing a confidential and a synthetic (as post synth object) dataset
generate_var <- function(p, corr, n, prop = 1){
  # relationship reversed (between V1 and everything else)
  vcov = matrix(corr, p, p) + diag(1-corr, p, p) # vcov matrix, p on all off-diagonal variables
  means = rep(0,p) # mean vector
  conf4 = as_tibble(mvrnorm(means, vcov, n = n), n=p)
  vcov_new = vcov
  vcov_new[1,2:p] = -rep(corr, p-1)
  vcov_new[2:p,1] = -rep(corr, p-1)
  if (prop != 1){
    synth4 = generate_postsynth(rbind(sample_n(as_tibble(mvrnorm(means, vcov_new, n = n), n=p), size=n*prop, replace = FALSE),
                                      as_tibble(mvrnorm(means, vcov, n = n*(1-prop)), n=p)))
  } else{
    synth4 = generate_postsynth(rbind(sample_n(as_tibble(mvrnorm(means, vcov_new, n = n), n=p), size=n*prop, replace = FALSE)))
  }
  
  list(conf4, synth4)
}

# function returns the discriminator object or discriminant-based metric values with 
# propensity scores for synthetic and confidential data
# synth: synthetic dataset as a postsynth object
# conf: confidential dataset
# mod: model used to fit discriminant-based metrics
# type: model type label for visualization
# name: describes synthetic dataset
# param: parameter to tune, if any
# return: if "Standard," returns metric values. Else, return discriminator.
make_disc <- function(synth, conf, mod, type, name, param = "none", return = "Standard") {
  rpart_rec <- recipe(.source_label ~ ., data = discrimination(synth, conf)$combined_data)
  
  if (type == "lr"){ # for logistic regression
    rpart_rec <- rpart_rec %>%
      step_interact(terms = ~ all_numeric_predictors():all_numeric_predictors()) %>% # add first order interaction terms
      step_normalize(all_predictors()) 
  }
  if (param == "none"){ # we aren't doing any tuning
    # set up discriminator
    d <- discrimination(synth, conf) %>%
      add_propensities(
        recipe = rpart_rec,
        spec = mod
      ) %>%
      add_discriminator_auc() %>%
      add_specks() %>%
      add_pmse() %>%
      add_pmse_ratio(times = 25)
  }
  else{ # we are tuning
    if (param == "penalty"){ # lasso model
      rec <- rpart_rec %>%
        step_interact(terms = ~ all_numeric_predictors():all_numeric_predictors()) %>% # add first order interaction terms
        step_normalize(all_predictors()) 
      grid <- grid_regular(penalty(), levels = 10)
    }
    else if (param == "cp"){ # tree model
      rec <- rpart_rec 
      grid <- grid_regular(cost_complexity(), levels= 10)
    }
    # set up discriminator
    d <- discrimination(synth, conf) %>%
      add_propensities_tuned(
        grid = grid, 
        recipe = rec,
        spec = mod
      ) %>%
      add_discriminator_auc() %>%
      add_specks() %>%
      add_pmse() %>%
      add_pmse_ratio(times = 25)
  }
  if (return == "Standard"){
    # extract the metrics
    return(c(d$pmse$.pmse[1], d$specks$.specks[1], d$discriminator_auc$.estimate[1], d$pmse$.pmse_ratio[1],#train
             d$pmse$.pmse[2], d$specks$.specks[2], d$discriminator_auc$.estimate[2], d$pmse$.pmse_ratio[2],#test
             type, name))
  } else{
    return(d)
  }
}


################################ CASE STUDY ################################
# function processes and returns input dataframe
#   df: inptu dataframe
#   var: "race" or "pov"
#   return: processed dataframe
process_data <- function(df, var){
  df <- as.data.frame(unclass(df),stringsAsFactors=TRUE)
  df$pov <- as.factor(df$pov)
  df$black <- as.factor(df$black)
  if (var == "race"){
    order = c("black", "pov", "SEX", "EMPSTAT", "EDUC", "HCOVANY", "VETSTAT", "FTOTINC", "AGE")
  }
  else if (var == "pov"){
    order = c("pov", "EMPSTAT", "EDUC", "black", "SEX", "HCOVANY", "VETSTAT", "FTOTINC", "AGE")
  }
  else{
    stop("var required to be one of 'race'', 'pov'")
  }
  df <- df %>%
    dplyr::select(order)
  df
}

# function to succinctly synthesize data in a specified synthesis order
#   start_data: starting data as provided to visit_sequence() and roadmap()
#   conf_data: confidential data as provided to roadmap()
#   var_order: synthesis order as provided to visit_sequence()
#   return: synthesized dataset as a postsynth object
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

# function to permute a given proportion of entries in each column of a datatable
#   dt: data table to permute
#   k: proportion of entries to permute
#   return: permuted data table
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

# function to synthesize datasets with different utilites by race
#   data_mi: the confidential data
#   var: the grouping variable
#   return: list of 5 synthetic datasets (postsynth objects):
#     synthesized from unchanged confidential dataset
#     synthesized from confidential dataset in which 25% of the entries for black people are permuted
#     synthesized from confidential dataset in which 75% of the entries for black people are permuted
#     synthesized from confidential dataset in which 25% of the entries for White people are permuted
#     synthesized from confidential dataset in which 75% of the entries for White people are permuted
make_all_synth_datasets <- function(data_mi, var){
  if (var == "race"){
    # set synthesis order manually
    order = c("pov", "SEX", "EMPSTAT", "EDUC", "HCOVANY", "VETSTAT", "FTOTINC", "AGE")
    
    ########## synthesize (race, from all input data) ##########
    conf_data = data_mi
    start_data = data.frame("black" = data_mi[,"black"])
    synth_black_good <- synth(start_data, conf_data, order)
    
    ######### synthesize (race, permute 25% of the entries in each column of black individuals) ##########
    black_conf <- data_mi %>%
      filter(black == 1)
    sample_black_conf <- permute_entries(black_conf, 0.25)
    noblack_conf = data_mi %>%
      filter(black == 0)
    conf_data <- rbind(sample_black_conf, noblack_conf)
    start_data = data.frame("black" = as.factor(c(rep(0, nrow(noblack_conf)), 
                                                  rep(1, nrow(data_mi)-nrow(noblack_conf)))))
    synth_black_med <- synth(start_data, conf_data, order)
    
    ######### synthesize (race, permute 75% of the entries in each column black individuals) ##########
    black_conf <- data_mi %>%
      filter(black == 1)
    sample_black_conf <- permute_entries(black_conf, .75)
    noblack_conf = data_mi %>%
      filter(black == 0)
    conf_data <- rbind(sample_black_conf, noblack_conf)
    start_data = data.frame("black" = as.factor(c(rep(0, nrow(noblack_conf)), 
                                                  rep(1, nrow(data_mi)-nrow(noblack_conf)))))
    synth_black_bad <- synth(start_data, conf_data, order)
    
    ######### synthesize (race, permute 25% of the entries in each column of White individuals) ##########
    black_conf <- data_mi %>%
      filter(black == 1)
    noblack_conf = data_mi %>%
      filter(black == 0)
    sample_noblack_conf <- permute_entries(noblack_conf, 0.25)
    conf_data <- rbind(sample_noblack_conf, black_conf)
    start_data = data.frame("black" = as.factor(c(rep(0, nrow(noblack_conf)), 
                                                  rep(1, nrow(data_mi)-nrow(noblack_conf)))))
    synth_noblack_med <- synth(start_data, conf_data, order)
    
    ######### synthesize (race, permute 75% of the entries in each column of White individuals) ##########
    black_conf <- data_mi %>%
      filter(black == 1)
    noblack_conf = data_mi %>%
      filter(black == 0)
    sample_noblack_conf <- permute_entries(noblack_conf, .75)
    conf_data <- rbind(sample_noblack_conf, black_conf)
    start_data = data.frame("black" = as.factor(c(rep(0, nrow(noblack_conf)), 
                                                  rep(1, nrow(data_mi)-nrow(noblack_conf)))))
    synth_noblack_bad <- synth(start_data, conf_data, order)
    
    return(list(synth_black_good, synth_black_med, synth_black_bad, synth_noblack_med, synth_noblack_bad))
  } 
  else if (var == "pov"){
    # variable ordering: 
    dmy <- dummyVars(" ~ .", data = data_mi)
    trsf <- data.frame(predict(dmy, newdata = data_mi))
    cor <- data.frame(cor(trsf))
    rownames(cor)[order(abs(cor$pov.1), decreasing = TRUE)]
    
    # based on this (select variables most highly correlated with pov first), and doing categorical variables first, order is:
    order = c("EMPSTAT", "EDUC", "black", "SEX", "HCOVANY", "VETSTAT", "FTOTINC", "AGE")
    
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
    
    return(list(synth_pov_good, synth_pov_med, synth_pov_bad, synth_nopov_med, synth_nopov_bad))
  } 
  else{
    stop("var required to be one of 'race'', 'pov'")
  }
}

# function to calculate discriminant-based metrics using a decision tree
#   discriminator: the input discriminator object
#   return: the discriminator object with propensities scores
calc_disc <- function(discriminator){
  # Evaluate discriminant-based metrics on the data using tree-based model
  tree_mod <- decision_tree(cost_complexity = tune()) %>%
    set_mode(mode = "classification") %>%
    set_engine(engine = "rpart")
  
  rpart_rec <- recipe(.source_label ~ ., data = discriminator$combined_data)
  
  grid = grid_regular(cost_complexity(), levels= 10)
  
  # set up discriminator
  d <- discriminator %>%
    add_propensities_tuned(
      grid = grid, 
      recipe = rpart_rec,
      spec = tree_mod
    )
  d
}

# function to fit two-model approach (create separate models for majority and minority group)
#   synth: the synthetic dataset
#   conf: the confidential dataset
#   calc_disc: function to fit the discriminator
#   synth_type: label to be included in output
#   var: the grouping variable
#   return: dataframe of AUC, SPECKS, synth_type, and group
fit_two_model <- function(synth, conf, calc_disc, synth_type, var){
  if (var == "race"){
    # black disc
    synth_black = synth[synth$black == 1,]
    conf_black = conf[conf$black == 1,]
    disc = discrimination(synth_black, conf_black)
    res = calc_disc(disc)
    
    res <- res %>%
      add_discriminator_auc() %>%
      add_specks()
    
    auc_black = res$discriminator_auc$.estimate[2]
    specks_black = res$specks$.specks[2]
    
    # no black disc
    synth_noblack = synth[synth$black == 0,]
    conf_noblack = conf[conf$black == 0,]
    disc = discrimination(synth_noblack, conf_noblack)
    res = calc_disc(disc)
    
    res <- res %>%
      add_discriminator_auc() %>%
      add_specks()
    
    auc_noblack = res$discriminator_auc$.estimate[2]
    specks_noblack = res$specks$.specks[2]
    
    return(data.frame(aucs = c(auc_black, auc_noblack), 
                      speckss = c(specks_black, specks_noblack), 
                      synth_type = rep(synth_type, 2),
                      black_class = c("black", "noblack")))
  } 
  else if (var == "pov"){
    # pov disc
    synth_pov = synth[synth$pov == 1,]
    conf_pov = conf[conf$pov == 1,]
    disc = discrimination(synth_pov, conf_pov)
    res = calc_disc(disc)
    
    res <- res %>%
      add_discriminator_auc() %>%
      add_specks()
    
    auc_pov = res$discriminator_auc$.estimate[2]
    specks_pov = res$specks$.specks[2]
    
    # no pov disc
    synth_nopov = synth[synth$pov == 0,]
    conf_nopov = conf[conf$pov == 0,]
    disc = discrimination(synth_nopov, conf_nopov)
    res = calc_disc(disc)
    
    res <- res %>%
      add_discriminator_auc() %>%
      add_specks()
    
    auc_nopov = res$discriminator_auc$.estimate[2]
    specks_nopov = res$specks$.specks[2]
    
    return(data.frame(aucs = c(auc_pov, auc_nopov), 
                      speckss = c(specks_pov, specks_nopov), 
                      synth_type = rep(synth_type, 2),
                      pov_class = c("pov", "nopov")))
  } else{
    stop("var required to be one of 'race'', 'pov'")
  }
  
}

# Given a discriminator, get metrics for one- and two-model approach
#   synth: the synthetic dataset
#   data_mi: the confidential dataset
#   calc_disc: function to fit the discriminator
#   synth_type: label to be included in output
#   var: the grouping variable
#   return: list containing dataframes of AUC, SPECKS, synth_type, and group for one- and two-model approaches
get_nums <- function(synth, data_mi, calc_disc, synth_type, var){
  if (var == "race"){
    disc = discrimination(synth, data_mi)
    res = calc_disc(disc)
    res <- res %>%
      add_discriminator_auc(group = c("black")) %>%
      add_specks(group = c("black"))
    auc <- res$discriminator_auc$.estimate[3:4]
    specks <- res$specks$.specks[3:4]
    one_model_df <- data.frame(aucs = auc, 
                               speckss = specks, 
                               synth_type = rep(synth_type, 2),
                               black_class = c("noblack", "black"))
    
    two_model_df <- fit_two_model(synth, data_mi, calc_disc, synth_type, var)
    
    return(list(one_model_df, two_model_df))
  } else if (var == "pov"){
    disc = discrimination(synth, data_mi)
    res = calc_disc(disc)
    res <- res %>%
      add_discriminator_auc(group = c("pov")) %>%
      add_specks(group = c("pov"))
    auc <- res$discriminator_auc$.estimate[3:4]
    specks <- res$specks$.specks[3:4]
    one_model_df <- data.frame(aucs = auc, 
                               speckss = specks, 
                               synth_type = rep(synth_type, 2),
                               pov_class = c("nopov", "pov"))
    
    two_model_df <- fit_two_model(synth, data_mi, calc_disc, synth_type, var)
    
    return(list(one_model_df, two_model_df))
  } else{
    stop("var required to be one of 'race'', 'pov'")
  }
}

# function to fit the pMSE-ratio with the denominator correction. 
#   synth: the synthetic dataset
#   data_mi: the confidential dataset
#   calc_disc: function used to calculate the discriminator
#   times: the number of bootstrap samples to draw to contribute to the denominator
#   minority_size: the number of samples in the minority group
corrected_pmse_ratio <- function(synth, data_mi, calc_disc, times, minority_size, var){
  if (var == "race"){
    # create no black sub-dataset
    conf_noblack <- data_mi %>%
      filter(black == 0)
    synth_noblack <- synth %>%
      filter(black == 0)
    disc_noblack <- discrimination(synth_noblack, conf_noblack)
    
    # fit model wtih pmse to our no black sub-dataset
    res_noblack <- calc_disc(disc_noblack)
    # fit pmse ratio using size of minority group
    obj <- res_noblack %>%
      add_pmse() %>%
      add_pmse_ratio_scaled(size = minority_size, times = times)
    res_noblack <- obj[[1]]
    denoms_noblack <- obj[[2]]
    
    # create black sub-dataset
    conf_black <- data_mi %>%
      filter(black == 1)
    synth_black <- synth %>%
      filter(black == 1)
    disc_black <- discrimination(synth_black, conf_black)
    
    # fit model with pmse to our no black sub-dataset
    res_black <- calc_disc(disc_black)
    
    # fit pmse ratio using size of minority group
    obj <- res_black %>%
      add_pmse() %>%
      add_pmse_ratio_scaled(size = minority_size, times = times)
    res_black <- obj[[1]]
    denoms_black <- obj[[2]]
    
    return(list(res_noblack$pmse$.pmse_ratio[2], res_black$pmse$.pmse_ratio[2], 
                res_noblack$pmse$.pmse[2]/denoms_noblack,
                res_black$pmse$.pmse[2]/denoms_black))
  } else if (var == "pov"){
    # create no poverty sub-dataset
    conf_nopov <- data_mi %>%
      filter(pov == 0)
    synth_nopov <- synth %>%
      filter(pov == 0)
    disc_nopov <- discrimination(synth_nopov, conf_nopov)
    
    # fit model wtih pmse to our no poverty sub-dataset
    res_nopov <- calc_disc(disc_nopov)
    # fit pmse ratio using size of minority group
    obj <- res_nopov %>%
      add_pmse() %>%
      add_pmse_ratio_scaled(size = minority_size, times = times)
    res_nopov <- obj[[1]]
    denoms_nopov <- obj[[2]]
    
    # create poverty sub-dataset
    conf_pov <- data_mi %>%
      filter(pov == 1)
    synth_pov <- synth %>%
      filter(pov == 1)
    disc_pov <- discrimination(synth_pov, conf_pov)
    
    # fit model wtih pmse to our no poverty sub-dataset
    res_pov <- calc_disc(disc_pov)
    
    # fit pmse ratio using size of minority group
    obj <- res_pov %>%
      add_pmse() %>%
      add_pmse_ratio_scaled(size = minority_size, times = times)
    res_pov <- obj[[1]]
    denoms_pov <- obj[[2]]
    
    return(list(res_nopov$pmse$.pmse_ratio[2], res_pov$pmse$.pmse_ratio[2], 
                sd(res_nopov$pmse$.pmse[2]/denoms_nopov)/sqrt(length(denoms_nopov)-1),
                sd(res_pov$pmse$.pmse[2]/denoms_pov)/sqrt(length(denoms_pov)-1)))
  } else{
    stop("var required to be one of 'race'', 'pov'")
  }
}

# function to fit the pMSE-ratio without the denominator correction. 
#   synth: the synthetic dataset
#   data_mi: the confidential dataset
#   calc_disc: function used to calculate the discriminator
#   times: the number of bootstrap samples to draw to contribute to the denominator
pmse_ratio <- function(synth, data_mi, calc_disc, times, var){
  if (var == "race"){
    # create no black sub-dataset
    conf_noblack <- data_mi %>%
      filter(black == 0)
    synth_noblack <- synth %>%
      filter(black == 0)
    disc_noblack <- discrimination(synth_noblack, conf_noblack)
    
    # fit model wtih pmse to our no black sub-dataset
    res_noblack <- calc_disc(disc_noblack)
    # fit pmse ratio using size of minority group
    obj <- res_noblack %>%
      add_pmse() %>%
      add_pmse_ratio(times = times)
    res_noblack <- obj$pmse$.pmse_ratio[2]
    
    # create black sub-dataset
    conf_black <- data_mi %>%
      filter(black == 1)
    synth_black <- synth %>%
      filter(black == 1)
    disc_black <- discrimination(synth_black, conf_black)
    
    # fit model with pmse to our no black sub-dataset
    res_black <- calc_disc(disc_black)
    
    # fit pmse ratio using size of minority group
    obj <- res_black %>%
      add_pmse() %>%
      add_pmse_ratio(times = times)
    res_black <- obj$pmse$.pmse_ratio[2]
    
    return(c(res_noblack, res_black))
  } else if (var == "pov"){
    # create no poverty sub-dataset
    conf_nopov <- data_mi %>%
      filter(pov == 0)
    synth_nopov <- synth %>%
      filter(pov == 0)
    disc_nopov <- discrimination(synth_nopov, conf_nopov)
    
    # fit model wtih pmse to our no poverty sub-dataset
    res_nopov <- calc_disc(disc_nopov)
    # fit pmse ratio using size of minority group
    obj <- res_nopov %>%
      add_pmse() %>%
      add_pmse_ratio(times = times)
    res_nopov <- obj$pmse$.pmse_ratio[2]
    
    # create poverty sub-dataset
    conf_pov <- data_mi %>%
      filter(pov == 1)
    synth_pov <- synth %>%
      filter(pov == 1)
    disc_pov <- discrimination(synth_pov, conf_pov)
    
    # fit model wtih pmse to our no poverty sub-dataset
    res_pov <- calc_disc(disc_pov)
    
    # fit pmse ratio using size of minority group
    obj <- res_pov %>%
      add_pmse() %>%
      add_pmse_ratio(times = times)
    res_pov <- obj$pmse$.pmse_ratio[2]
    
    return(c(res_nopov, res_pov))
  } else{
    stop("var required to be one of 'race'', 'pov'")
  }
}

