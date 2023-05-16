# ######################################################################
# title model_helper_functions.R
# author  Lydia Haile
# purpose helper functions for model runs
########################################################################


run_malaria_model<- function(output, folder, stochastic_run){
  
  #' Prep inputs (without parameter updates)-- not explicitly passing site_data object in but this should be updated
  #' 
  #' @param output         list with site name, urban/rural grouping, and parameters to pass into cluster
  #'                       generated from prep_inputs function
  #' @param stochastic_run stochastic or central burden estimate run? 
  #'                       If this is a model run for stochastic burden estimates, adds an identifying column for the draw number corresponding to the posterior MCMC distribution.
  #' output: model outputs for site with provided parameters (saved into pre-specified folder)
  
  # run the model
  message('running the model')
  require(data.table)
  
  model<- malariasimulation::run_simulation(timesteps = output$param_list$timesteps,
                                            parameters = output$param_list) 
  
  model<- data.table(model)
  model[, site_name:= output$site_name]
  model[, urban_rural:=output$ur]
  model[, iso:= output$iso]
  
  if (stochastic_run== T){
  
 model[, draw:= output$stochastic_draw_number]
    
  }
  # save model runs somewhere
  message('saving the model')
  saveRDS(model, file= paste0(folder, 'raw_model_output_', output$site_name, '_', output$ur, '.RDS'))
}


