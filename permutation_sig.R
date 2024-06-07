#' Calculates permutation test-based confidence intervals
#'
#' @param conf_data the confidential input dataset -- column to permute will be taken from here
#' @param synth_data the synthetic dataset
#' @param col the column to permute
#' @param func function to apply (of form add_* from syntheval)
#' @param N the number of permutations to create
#' @param calc_disc function to run discriminator model
#' 
#' @return An *xN matrix where rows correspond to the test outputs of func
#' and each column corresponds to the output from one permutation
#' 
#' @export
#'
permutation_sig <- function(conf_data, synth_data, col, func, N, calc_disc, name) {
  colnames = paste0(col, levels(conf_data[,col]))
  result = data.frame(matrix(nrow = 0, ncol = length(colnames))) 
  colnames(result) = colnames
  
  for (i in 1:N){
    # Select randomly n rows of synthetic/confidential datasets, where n = number of people in poverty
    sample = permute(conf_data[,col])
    
    # Set up new pov indicator in confidential and synthetic data
    synth_sample = synth_data
    synth_sample[,col] = sample
    
    conf_sample = conf_data
    conf_sample[,col] = sample
    
    disc_sample = calc_disc(discrimination(synth_sample, conf_sample)) %>%
      func(group = c(col))
    
    # take last half of rows in dataframe
    toadd = data.frame(disc_sample[name])[(nrow(data.frame(disc_sample[name]))/2 + 1):nrow(data.frame(disc_sample[name])),
                                          ncol(data.frame(disc_sample[name]))]
    result[nrow(result)+1,] = toadd
  }
  return(result)
}