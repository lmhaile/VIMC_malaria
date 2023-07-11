# ######################################################################
# title model_helper_functions.R
# author  Lydia Haile
# purpose helper functions for model runs
########################################################################

#' run malaria model
#'
#' @param filepath filepath where model inputs are 

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
                              site_name,
                              '_',
                              ur,
                              '_',
                              iso,
                              '_', 
                              tag, 
                              '.RDS'))
  
  
}



run_malaria_model_rfp<- function(input){
  
  #' Prep inputs (without parameter updates)-- not explicitly passing site_data object in but this should be updated
  #' 
  #' @param input    list with parameter_filepath and output_folder
  #' @param tag tag to save input file with for identification
  #' output: model outputs for site with provided parameters (saved into pre-specified folder)
  
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
