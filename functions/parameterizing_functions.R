######################################################
## title    Parameterization functions
## author   Lydia Haile
## purpose  functions used to paramaterize VIMC model runs. Includes:
##            
######################################################

#' Set vaccine coverage using scene package
#'
#' @param   site          site data file
#' @param   change        would you like to set a certain coverage level? Boolean
#' @param   terminal_year year you would like to expand intervention coverage out to 
#' @param   rtss_target   data set with mortality inputs
#' @param   rtss_year     year for target
#' @returns site file with updated vaccine coverage values 

set_vaccine_coverage<- function(site, change, terminal_year, rtss_target, rtss_year){
  
  
  require(scene)
  
  group_var <- names(site$sites)
  
  # expand intervention years ----------------------------------------------------
  site$interventions <- site$interventions |>
    expand_interventions(max_year = terminal_year,
                         group_var = group_var)
  
  if (change== TRUE){
    site$interventions <- site$interventions |>
      set_change_point(sites = site$sites, 
                       var = "rtss_cov", 
                       year = rtss_year, 
                       target = rtss_target)
    
    # Linear scale up of coverage
    # if you would like coverage to scale up to a certain target
    site$interventions <- site$interventions |>
      linear_interpolate(vars = c("itn_use", "pmc_cov", "smc_cov", "rtss_cov"), 
                         group_var = group_var)
    
    
  }
  
  
  site$interventions <- site$interventions |>
    fill_extrapolate(group_var = group_var)
  
  
  site$interventions <- site$interventions |>
    add_future_net_dist(group_var = group_var)
  
  
  return(site)
}



#' Prep inputs for batch launch for central burden estimate GAVI runs
#'
#' @param site_data  dataset with site files for country
#' @param folder     folder to save input parameters
#' @param population population to run the model on
#' @param min_ages   minimum ages
#' @param max_ages   maximum ages
#' @returns list with site name, urban/rural grouping, iso code, and parameters to pass into cluster 


prep_inputs<- function(site_data, folder, population, min_ages, max_ages){
  
  
  # how many sites in this country?
  jobs<- nrow(site_data$sites)
  
  message(paste0('prepping ', jobs, ' jobs for model launch'))
  
  prep_site_data<- function(num){
    site<- site::single_site(site_file= site_data, index= num) 
    
    ## get site info
    site_name<- site$sites$name_1
    ur<- site$sites$urban_rural
    iso<- site$sites$iso3c
    message(paste0('prepping inputs for site ', site_name, ' ', ur))
    
    # pull parameters for this site
    params<- site::site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir= site$eir$eir[1],
      overrides = list(human_population= population),
      burnin = 5* 365
    )
    
    year<- 365
    # Set clinical incidence rendering
    params$clinical_incidence_rendering_min_ages = min_ages 
    params$clinical_incidence_rendering_max_ages = max_ages
    
    # Set severe incidence rendering
    params$severe_incidence_rendering_min_ages = min_ages 
    params$severe_incidence_rendering_max_ages = max_ages 
    
    # Set clinical incidence rendering
    params$clinical_incidence_rendering_min_ages =  min_ages
    params$clinical_incidence_rendering_max_ages = max_ages 
    
    # Set age group rendering
    params$age_group_rendering_min_ages = min_ages 
    params$age_group_rendering_max_ages = max_ages 
    
    inputs<- list('param_list'= params, 'site_name'= site_name, 'ur'= ur, 'iso'= iso)
    
    write_rds(inputs, paste0(folder, 'input_parameters/', site_name, '_', ur, '_', iso, '.rds'))
    
    print('saved input')
    
    return(paste0(site_name, '_', ur, '_', iso))
  }
  output<- lapply(c(1:jobs), prep_site_data)
}

#' Prep inputs for batch launch for central burden estimate GAVI runs
#'
#' @param site_data  dataset with site files for country
#' @param mort_dt    dataset with mortality inputs
#' @param folder     folder to save input parameters
#' @param population population to run the model on
#' @returns list with site name, urban/rural grouping, iso code, and parameters to pass into cluster


prep_inputs_rfp<- function(site_data, mort_dt, death_rate_matrix, folder){
  
  
  # how many sites in this country?
  jobs<- nrow(site_data$sites)
  
  message(paste0('prepping ', jobs, ' jobs for model launch'))
  
  prep_site_data<- function(num){
    site<- site::single_site(site_file= site_data, index= num) 
    
    ## get site info
    site_name<- site$sites$name_1
    ur<- site$sites$urban_rural
    iso<- site$sites$iso3c
    message(paste0('prepping inputs for site ', site_name, ' ', ur))
    
    # pull parameters for this site
    params<- site::site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir= site$eir$eir[1],
      overrides = list(human_population= 5000),
      burnin = 5* 365
    )
    
    year<- 365
    # Set clinical incidence rendering
    params$clinical_incidence_rendering_min_ages = c(seq(0, 99, by = 1))*year
    params$clinical_incidence_rendering_max_ages = c(seq(1, 99, by = 1)*year -1, 36500)
    
    # Set severe incidence rendering
    params$severe_incidence_rendering_min_ages = c(seq(0, 99, by = 1))*year
    params$severe_incidence_rendering_max_ages = c(seq(1, 99, by = 1)*year -1, 36500)
    
    # Set clinical incidence rendering
    params$clinical_incidence_rendering_min_ages = c(seq(0, 99, by = 1))*year
    params$clinical_incidence_rendering_max_ages = c(seq(1, 99, by = 1)*year -1, 36500)
    
    # Set age group rendering
    params$age_group_rendering_min_ages = c(seq(0, 99, by = 1))*year
    
    inputs<- list('param_list'= params, 'site_name'= site_name, 'ur'= ur, 'iso'= iso)
    
    write_rds(inputs, paste0(folder, 'input_parameters/', site_name, '_', ur, '_', iso, '.rds'))
    
    return(paste0(site_name, '_', ur, '_', iso))
  }
  output<- lapply(c(1:jobs), prep_site_data)
}


#' Pull stochastic parameters and calibrates for stochastic model runs
#' 
#' @param site List created by prep_inputs. 
#'                  List contains site name, urban/rural grouping, iso code, and parameters to pass into cluster
#' @param draws     The number of stochastic parameter draws you would like to pull
#' @returns list with site name, urban/rural grouping, iso code, and additional 'stoch_draws' list with stochastic draws equal to the number of draws requested.

prep_stochastic_inputs<- function(site, draws){
  
  
  format_stochastic_inputs<- function(x){
    
    message(paste0('pulling inputs for draw number ', x))
    
    param_draw<- site$param_list |>
      set_parameter_draw(x) |>
      set_equilibrium(init_EIR= site$param_list$init_EIR)
    
    stoch_inputs<- list('param_list'= param_draw,
                        'site_name' = site$site_name,
                        'ur' = site$ur,
                        'iso' = site$iso,
                        'stochastic_draw_number' = x)
    
    return(stoch_inputs)
  }
  
  output<- lapply(1:draws, format_stochastic_inputs)
  
  return(output)
}

#' Quick calibration for VIMC models using calibrate
#' 
#' @param site_data site file
#' @returns calibrated EIR estimate


quick_calibration<- function(site_data){
  
  require(cali)
  
  prev<- data.table(site_data$prevalence)
  target_pfpr <- prev[year == 2007, pfpr]
  
  # Define our simulation parameters, we need to add $timesteps for the calibration
  p <- malariasimulation::get_parameters(
    overrides = list(
      human_population = 5000,
      individual_mosquitoes = FALSE
    )
  )
  p$timesteps <- 7 * 365
  
  # Write our summary function. To match the target we need to
  # output average PfPr in year 3
  get_pfpr_year_3 <- function(simulation_output){
    year3 <- simulation_output[simulation_output$timestep > 6 * 365,]
    pfpr <- year3$n_detect_730_3650 / year3$n_730_3650
    average_pfpr <- mean(pfpr)
    return(average_pfpr)
  }
  
  # Run a test simulation to check our target function does what it should!
  simulation <- malariasimulation::run_simulation(
    timesteps = p$timesteps,
    parameters = p
  )
  # Does this look right? Does it match the type of target data?
  get_pfpr_year_3(simulation)
  
  # When we are happy we can run the calibration
  set.seed(1234)
  calibration <- cali::calibrate(target = target_pfpr,
                                 summary_function = get_pfpr_year_3,
                                 parameters = p,
                                 tolerance = 0.005,
                                 low = 0.02, high = 20)
  
  message("Calibrated EIR estimate: ", calibration)
  return(calibration)
  
}

