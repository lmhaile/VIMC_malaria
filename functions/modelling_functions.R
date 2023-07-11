# ######################################################################
# title model_helper_functions.R
# author  Lydia Haile
# purpose helper functions for model runs
########################################################################

#' Run malariasimulation model for VIMC 
#'
#' @param filepath filepath where model inputs are (list with model parameters and demographic info, produced by prep_model_launch)
#' @returns model output saved in pre-specified VIMC_files folder
#' 

run_malaria_model<- function(filepath) {
  
  model_input<- readRDS(filepath)
  
  params<-model_input$param_list
  params$progress_bar<- TRUE
  timesteps<<- model_input$param_list$timesteps
  

  scenario<- model_input$scenario
  tag<- model_input$tag
  
  
  model<- malariasimulation::run_simulation(timesteps = params$timesteps,
                                            parameters = params) 
  
  model<- data.table(model)
  model[, site_name:= model_input$site_name]
  model[, urban_rural:= model_input$ur]
  model[, iso:= model_input$iso]
  
  
  # save model runs somewhere
  message('saving the model')
  saveRDS(model, file= paste0('Q:/VIMC_files/central_estimates/', 
                              scenario, 
                              'raw_model_output/raw_model_output_', 
                              model_input$site_name,
                              '_',
                              model_input$ur,
                              '_',
                              model_input$iso,
                              '_', 
                              tag, 
                              '.RDS'))
  
  
}


#' Prep inputs (without parameter updates)-- not explicitly passing site_data object in but this should be updated
#' 
#' @param input    list with parameter_filepath and output_folder
#' @param tagging tag to save input file with for identification
#' @returns model outputs for site with provided parameters (saved into pre-specified folder) 

run_malaria_model_rfp<- function(input){
  
  
  # run the model
  message('running the model')
  
  
  # read in model inputs
  params<- readRDS(input$parameter_filepath)
  tagging<- input$tag
  
  print(input$parameter_filepath)
  print(tagging)
  
  if(is.list(params)== TRUE){
    
    message('succesfully loaded input parameters')
    
  }
  message('running the model')
  
  model<- malariasimulation::run_simulation(timesteps = 365*65,
                                            parameters = params) 
  
  message('finished running the model')
  
  model<- data.table(model)
  fold<- input$output_folder
  
  # save model runs somewhere
  message('saving the model')
  saveRDS(model, file= paste0(fold, 'raw_model_output_rfp_', tagging, '.rds'))
  
  message('successfully saved the model')
}
