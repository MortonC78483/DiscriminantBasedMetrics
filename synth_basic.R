#' Basic synthesis of data from variable ordering passed in
#'
#' @param start_data starting data for synthesis
#' @param conf_data confidential data for synthesis
#' @param var_order a vector of variable names dictating the visit order
#' (does not include variables in start_data)
#'
#' @return postsynth object from starting and confidential data
#' 
#' @family Synthesis
#' 
#' @export

# function that takes start data and confidential data, and synthesizes dataset
synth <- function(start_data, conf_data, var_order){
  visit_sequence = visit_sequence(var_order, start_data, type = "manual")
  roadmap = roadmap(conf_data, start_data, visit_sequence)
  recipe = construct_recipes(roadmap = roadmap)
  
  # use basic tree models, for either classification or regression, as appropriate
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