################################################################################
## title    Parameterization functions
## author   Lydia Haile
## purpose  functions used to paramaterize VIMC model runs. Includes:
################################################################################

#' Set vaccine coverage using scene package
#'
#' @param   site          site data file
#' @param   scenario      'baseline' or 'intervention'
#' @param   terminal_year year you would like to expand intervention coverage out to 
#' @param   rtss_target   data set with mortality inputs
#' @param   rtss_year     year for target
#' @returns site file with updated vaccine coverage values 
set_vaccine_coverage<- function(site, 
                                scenario, 
                                terminal_year, 
                                rtss_target,
                                rtss_year){
  
  require(scene)
  group_var <- names(site$sites)
  
  # expand intervention years --------------------------------------------------
  site$interventions <- site$interventions |>
    expand_interventions(max_year = terminal_year,
                         group_var = group_var)
  
  if (scenario== 'intervention'){
    
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


set_vimc_mortality<- function(site){
  # bind on a year for youngest age group
  youngest<- data.table(site$demography)[age_upper== min(site$demography$age_upper)]
  mort_dt<- rbind(youngest, mort)
  site$demography<- mort_dt
  
  return(site)
}
#' Format input parameters for VIMC sites
#'
#' @param   site               site file for one site
#' @param   mortality          if TRUE, incorporate VIMC mortality inputs
#' @param   population         human population  value
#' @param   scenario           'baseline' or 'intervention'
#' @param   min_ages           minimum ages (in timesteps) for incidence/ prevalence outputs
#' @param   max_ages           maximum ages (in timesteps) for incidence/ prevalence outputs
#' @param   description                description of model runs you would like to carry out
#' @inheritParams set_vaccine_coverage
#' @returns list with input parameters and demographic info, as inputs for run_malaria_model
prep_site<- function(site,
                     save_dir= 'Q:/VIMC/central_estimates/input_parameters/',
                     vimc_mortality, 
                     population, 
                     scenario, 
                     min_ages, 
                     max_ages, 
                     description){
  
  # get site info
  site_name<- site$sites$name_1
  ur<- site$sites$urban_rural
  iso<- site$sites$iso3c
  
  
  # create directories to save outputs  -------
  if(dir.exists(paste0(save_dir, iso))== FALSE){
    dir.create(paste0(save_dir, iso))
  }
  
  if(dir.exists(paste0(save_dir, iso, '/', description))== FALSE){
    dir.create(paste0(save_dir, iso, '/', description))
  }
  
  message(paste0('prepping inputs for site ', site_name, ' ', ur))
  
  if(vimc_mortality== TRUE){
  site<- set_vimc_mortality(site)
  }

  if (scenario== 'baseline'){
    # # expand scenario out to 2050
    site <- set_vaccine_coverage(
      site,
      scenario = scenario,
      terminal_year = 2050,
      rtss_target = NULL,
      rtss_year = NULL
    )
  } else {
    site <- set_vaccine_coverage(
      site,
      scenario = scenario,
      terminal_year = 2050,
      rtss_target = 0.8,
      rtss_year = 2023
    )
  }
  # pull parameters for this site
  params<- site::site_parameters(
    interventions = site$interventions,
    demography = site$demography,
    vectors = site$vectors,
    seasonality = site$seasonality,
    eir= site$eir$eir[1],
    overrides = list(human_population= population)
  )
  
  year<- 365
  
  # Set clinical incidence rendering
  params$clinical_incidence_rendering_min_ages = min_ages 
  params$clinical_incidence_rendering_max_ages = max_ages
  
  # Set severe incidence rendering
  params$severe_incidence_rendering_min_ages = min_ages 
  params$severe_incidence_rendering_max_ages = max_ages 
  
  # Set clinical incidence rendering
  params$clinical_incidence_rendering_min_ages = min_ages
  params$clinical_incidence_rendering_max_ages = max_ages 
  
  # Set age group rendering
  params$age_group_rendering_min_ages = min_ages 
  params$age_group_rendering_max_ages = max_ages 
  
  inputs<- list('param_list'= params, 
                'site_name'= site_name, 
                'ur'= ur, 
                'iso'= iso,
                'scenario'= scenario,
                'description'= description)
  
  write_rds(
    inputs,
    paste0('Q:/VIMC/central_estimates/input_parameters/',
           iso, '/', 
           description, '/', 
           site_name, '_', 
           ur, '_', 
           scenario,
           '.rds'
    )
  )
  
}

#' Format input parameters for VIMC sites
#'
#' @param   site_data          site data file for a particular country
#' @param   mortality          if TRUE, incorporate VIMC mortality inputs
#' @param   scenario           'baseline' or 'intervention'
#' @param   min_ages           minimum ages (in timesteps) for incidence/ prevalence outputs
#' @param   max_ages           maximum ages (in timesteps) for incidence/ prevalence outputs
#' @param   description                description of model runs you would like to carry out
#' @inheritParams set_vaccine_coverage
#' @returns list with input parameters and demographic info, as inputs for run_malaria_model
prep_country<- function(iso, 
                        save_dir= 'Q:/VIMC/central_estimates/input_parameters/',
                        vimc_mortality,
                        population,
                        scenario, 
                        min_ages, 
                        max_ages, 
                        description){
  
  require(foresite)
  
  site_data<- get(iso)
  
  # create directories to save outputs  -------
  if(dir.exists(paste0(save_dir, iso))== FALSE){
    dir.create(paste0(save_dir, iso))
  }
  
  if(dir.exists(paste0(save_dir, iso, '/', description))== FALSE){
    dir.create(paste0(save_dir, iso, '/', description))
  }
  
  sites<- nrow(site_data$sites)
  
  for(num in c(1:sites)){
    message(paste0('prepping site ', num))
    
    site<- site::single_site(site_file= site_data, index= num) 
    prep_site(site, 
              vimc_mortality= vimc_mortality, 
              scenario= scenario, 
              population= population, 
              min_ages= min_ages, 
              max_ages= max_ages, 
              description= description)
    
  }
  
return(message(paste0('prepped outputs for model run: ', description, ', country: ', iso)))
  
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

