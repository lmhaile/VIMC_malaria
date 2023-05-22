# ######################################################################
# title model_helper_functions.R
# author  Lydia Haile
# purpose helper functions for model runs
########################################################################


run_malaria_model<- function(identifier, folder, stochastic_run){
  
  #' Prep inputs (without parameter updates)-- not explicitly passing site_data object in but this should be updated
  #' 
  #' @param identifier     site name, urban/rural grouping, and parameters to pass into cluster
  #'                       generated from prep_inputs function
  #' @param stochastic_run stochastic or central burden estimate run? 
  #'                       If this is a model run for stochastic burden estimates, adds an identifying column for the draw number corresponding to the posterior MCMC distribution.
  #'                       
  #' @param folder        folder to save model output. 
  #' output: model outputs for site with provided parameters (saved into pre-specified folder)
  
  # run the model
  message('running the model')
  require(data.table)
  
  # read in model inputs
  input<- readRDS(paste0(folder,'input_parameters/', identifier, '.rds'))
  
  model<- malariasimulation::run_simulation(timesteps = input$param_list$timesteps,
                                            parameters = input$param_list) 
  
  model<- data.table(model)
  model[, site_name:= input$site_name]
  model[, urban_rural:=input$ur]
  model[, iso:= input$iso]
  
  if (stochastic_run== T){
  
 model[, draw:= input$stochastic_draw_number]
    
  }
  
#  model <- model |> 
#    select(timestep, 
#           ft,
#           iso, 
#           site_name,
#           urban_rural,
#           run,
#           contains("n_inc_clin"),   # clinical incidence
#           contains("n_inc_severe"), # severe incidence
#           contains("n_age")         # population
#    )
  
  # save model runs somewhere
  message('saving the model')
  saveRDS(model, file= paste0(folder, 'raw_model_output/raw_model_output_', input$site_name, '_', input$ur, '.RDS'))
}

random_walk <- function(x, n_steps) {
  ret <- numeric(n_steps)
  for (i in seq_len(n_steps)) {
    x <- rnorm(1, x)
    ret[[i]] <- x
  }
  ret
}
