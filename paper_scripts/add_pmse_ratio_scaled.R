#' Add pMSE ratio to discrimination object
#'
#' @param discrimination A discrimination with added propensities
#' @param split A logical for if the metric should be calculated separately for 
#' the training/testing split. Defaults to TRUE.
#' @param prop The proportion of data to be retained for modeling/analysis in 
#' the training/testing split. The sampling is stratified by the original and
#' synthetic data.
#' @param size size of dataset to resample from in each bootstrap sample, default to normal pMSE ratio
#' @param times The number of bootstrap samples.
#'
#' @return A discrimination with pMSE and pMSE ratio, and the denominators used (for uncertainty quantification)
#' 
#' @family Utility metrics
#' 
#' @export
add_pmse_ratio_scaled <- function(discrimination, split = TRUE, prop = 3 / 4, size = NULL, times) {
  
  if (is.null(discrimination$pmse)) {
    
    stop("Error: discrimination must have a pmse. Use add_pmse() before add_pmse_ratio()")
    
  }
  
  calc_pmse <- function(propensities) {
    
    # calculate the expected propensity
    prop_synthetic <- propensities %>%
      dplyr::summarize(
        n_synthetic = sum(.data$.source_label == "synthetic"),
        n_total = dplyr::n()
      ) %>%
      dplyr::mutate(prop_synthetic = .data$n_synthetic / .data$n_total) %>%
      dplyr::pull("prop_synthetic")
    
    propensities_vec <- propensities %>%
      dplyr::pull(".pred_synthetic")
    
    # calculate the observed pMSE
    pmse <- mean((propensities_vec - prop_synthetic) ^ 2)
    
    return(pmse)
  }
  
  pmse_null_overall <- vector(mode = "numeric", length = times)
  pmse_null_training <- vector(mode = "numeric", length = times)
  pmse_null_testing <- vector(mode = "numeric", length = times)
  
  for (a in seq_along(pmse_null_overall)) {
    
    # original code if no size passed
    if (is.null(size)){
      # bootstrap sample original observations to equal the size of the combined 
      # data
      # append the original labels so the proportions match
      bootstrap_sample <- dplyr::bind_cols(
        discrimination$combined_data %>%
          dplyr::filter(.data$.source_label == "original") %>%
          dplyr::slice_sample(n = nrow(discrimination$combined_data), replace = TRUE) %>%
          dplyr::select(-".source_label"),
        discrimination$combined_data %>%
          dplyr::select(".source_label")
      )
    }
    # new code if size passed in
    else{
      # bootstrap sample original observations to equal the size of the combined 
      # data
      # append the original labels so the proportions match
      combined_sample <- discrimination$combined_data %>%
        dplyr::filter(.data$.source_label == "original") %>%
        sample_n(size, replace = FALSE)
      
      bootstrap_sample <- combined_sample %>%
        dplyr::slice_sample(n = 2*size, replace = TRUE) %>%
        dplyr::select(-".source_label")
      bootstrap_sample$.source_label = c(rep("original", size), rep("synthetic", size))
    }
    
    if (split) {
      
      # make training/testing split
      data_split <- rsample::initial_split(
        data = bootstrap_sample,
        prop = prop,
        strata = ".source_label"
      )
      
      # fit the model from the pMSE on the bootstrap sample
      fitted_model <- parsnip::fit(
        discrimination$discriminator, 
        data = rsample::training(data_split)
      )
      
      # calculate the propensities
      propensities_df <- dplyr::bind_cols(
        stats::predict(fitted_model, new_data = discrimination$combined_data, type = "prob")[, ".pred_synthetic"],
        discrimination$combined_data
      ) %>%
        dplyr::mutate(
          .sample = dplyr::if_else(
            dplyr::row_number() %in% data_split$in_id, 
            true = "training", 
            false = "testing"
          )
        )
      
      # calculate the pmse for each bootstrap
      pmse_null_overall[a] <- calc_pmse(propensities_df)
      pmse_null_training[a] <- propensities_df %>%
        dplyr::filter(.data$.sample == "training") %>%
        calc_pmse()
      pmse_null_testing[a] <- propensities_df %>%
        dplyr::filter(.data$.sample == "testing") %>%
        calc_pmse()
      
    } else {
      
      # fit the model from the pMSE on the bootstrap sample
      fitted_model <- parsnip::fit(
        discrimination$discriminator, 
        data = bootstrap_sample
      )
      
      # calculate the propensities
      propensities_df <- dplyr::bind_cols(
        stats::predict(fitted_model, new_data = discrimination$combined_data, type = "prob")[, ".pred_synthetic"],
        discrimination$combined_data
      )
      
      # calculate the pmse for each bootstrap
      pmse_null_overall[a] <- calc_pmse(propensities_df)
      
    }
    
  }
  
  # find the mean of the bootstrapped pMSEs
  mean_null_pmse_overall <- mean(pmse_null_overall)
  mean_null_pmse_training <- mean(pmse_null_training)
  mean_null_pmse_testing <- mean(pmse_null_testing)
  
  # calculate the ratio for the training/testing split or overall data
  if (all(c("training", "testing") %in% discrimination$pmse$.source)) {
    
    pmse <- dplyr::bind_cols(
      discrimination$pmse,
      tibble::tibble(.null_pmse = c(mean_null_pmse_training, mean_null_pmse_testing))
    ) %>%
      dplyr::mutate(.pmse_ratio = .data$.pmse / .data$.null_pmse)
    denoms = pmse_null_testing
    
  } else {
    
    pmse <- dplyr::bind_cols(
      discrimination$pmse,
      tibble::tibble(.null_pmse = mean_null_pmse_overall)
    ) %>%
      dplyr::mutate(.pmse_ratio = .data$.pmse / .data$.null_pmse)
    denoms = pmse_null_training
  }
  
  discrimination$pmse <- pmse
  
  return(list(discrimination, denoms))
  
}