# ######################################################################
# title model_helper_functions.R
# author  Lydia Haile
# purpose helper functions for model runs
########################################################################


#' run_malaria_model<- function(identifier, folder, stochastic_run, tag){
#'   
#'   #' Prep inputs (without parameter updates)-- not explicitly passing site_data object in but this should be updated
#'   #' 
#'   #' @param identifier     site name, urban/rural grouping, and parameters to pass into cluster
#'   #'                       generated from prep_inputs function
#'   #' @param stochastic_run stochastic or central burden estimate run? 
#'   #'                       If this is a model run for stochastic burden estimates, adds an identifying column for the draw number corresponding to the posterior MCMC distribution.
#'   #'                       
#'   #' @param folder         folder to save model output. 
#'   #' @param tag            additional identifying clause to add to file name (if you are running several tests on the same site)
#'   #' output: model outputs for site with provided parameters (saved into pre-specified folder)
#'   
#'   # run the model
#'   message('running the model')
#'   require(data.table)
#'   
#'   # read in model inputs
#'   input<- readRDS(paste0(folder,'input_parameters/', identifier, '.rds'))
#'   
#'   params<-input$param_list
#'   params$progress_bar<- FALSE
#'   
#'   message('read inputs successfully')
#'   model<- malariasimulation::run_simulation(timesteps = params$timesteps,
#'                                             parameters = params) 
#'   
#'   model<- data.table(model)
#'   model[, site_name:= input$site_name]
#'   model[, urban_rural:=input$ur]
#'   model[, iso:= input$iso]
#'   
#'   if (stochastic_run== T){
#'   
#'  model[, draw:= input$stochastic_draw_number]
#'     
#'   }
#'   
#' #  model <- model |> 
#' #    select(timestep, 
#' #           ft,
#' #           iso, 
#' #           site_name,
#' #           urban_rural,
#' #           run,
#' #           contains("n_inc_clin"),   # clinical incidence
#' #           contains("n_inc_severe"), # severe incidence
#' #           contains("n_age")         # population
#' #    )
#'   
#'   # save model runs somewhere
#'   message('saving the model')
#'   saveRDS(model, file= paste0(folder, 'raw_model_output/raw_model_output_', identifier, '_', tag, '.RDS'))
#' }
#' 
#' 
#' 
#' run_malaria_model_rfp<- function(input){
#'   
#'   #' Prep inputs (without parameter updates)-- not explicitly passing site_data object in but this should be updated
#'   #' 
#'   #' @param input    list with parameter_filepath and output_folder
#'   #' @param tag tag to save input file with for identification
#'   #' output: model outputs for site with provided parameters (saved into pre-specified folder)
#'   
#'   # run the model
#'   message('running the model')
#'   
#' 
#'   # read in model inputs
#'   params<- readRDS(input$parameter_filepath)
#'   tagging<- input$tag
#'   
#'   print(input$parameter_filepath)
#'   print(tagging)
#'   
#'   if(is.list(params)== TRUE){
#'     message('succesfully loaded input parameters')
#'   }
#'   message('running the model')
#'   
#'   model<- malariasimulation::run_simulation(timesteps = 365*65,
#'                                             parameters = params) 
#'   
#'   message('finished running the model')
#'   
#'   model<- data.table(model)
#'   fold<- input$output_folder
#'   # save model runs somewhere
#'   message('saving the model')
#'   saveRDS(model, file= paste0(fold, 'raw_model_output_rfp_', tagging, '.rds'))
#'   
#'   message('successfully saved the model')
#' }



run_malaria_model<- function(model_input) {
  
  params<-model_input$param_list
  params$progress_bar<- TRUE
  timesteps<<- model_input$param_list$timesteps
  
  
  scenario<- model_input$scenario
  tag<- model_input$tag
  
  
  model<- malariasimulation::run_simulation(timesteps = params$timesteps,
                                            parameters = params) 
  
  model<- data.table(model)
  model[, site_name:= input$site_name]
  model[, urban_rural:=input$ur]
  model[, iso:= input$iso]
  
  
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

