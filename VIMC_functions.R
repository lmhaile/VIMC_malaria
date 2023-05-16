######################################################
## title    VIMC functions
## author   Lydia Haile
## purpose  functions used for VIMC malaria model runs
######################################################


# parameterizing inputs --------------------------------------------------------

set_vaccine_coverage<- function(site, change= TRUE, terminal_year, rtss_target, rtss_year){
  
  #' set vaccine coverage using scene package
  #'
  #' @param terminal_year year you would like to expand intervention coverage out to 
  #' @param change        would you like to set a certain coverage level? Boolean
  #' @param rtss_target   dataset with mortality inputs
  #' @param rtss_year     year for target
  #' @param site          site data file 
  #' 
  
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

prep_inputs<- function(site_data, mort_dt, death_rate_matrix){
  
  #' Prep inputs for batch launch for central burden estimate GAVI runs
  #'
  #' @param site_data dataset with site files for country
  #' @param mort_dt   dataset with mortality inputs
  #' output: list with site name, urban/rural grouping, iso code, and parameters to pass into cluster
  
  
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
      overrides = list(human_population= 10000) # what size population is appropriate?
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
    params$age_group_rendering_max_ages = c(seq(1, 99, by = 1)*year -1, 36500)
    
    # # set custom demography based on mortality inputs ------------------------
    # params<- set_demography(
    #   params,
    #   agegroups= unique(mort_dt$age_to),
    #   timesteps = unique(mort_dt$year)*365,
    #   deathrates = death_rate_matrix
    #   
    # )
    # 
    inputs<- list('param_list'= params, 'site_name'= site_name, 'ur'= ur, 'iso'= iso)
    return(inputs)
  }
  output<- lapply(c(1:jobs), prep_site_data)
}

prep_stochastic_inputs<- function(site, draws){
  
  #' Pull stochastic parameters and calibrates for stochastic model runs
  #'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
  #' @param site_data List created by prep_inputs. 
  #'                  List contains site name, urban/rural grouping, iso code, and parameters to pass into cluster
  #' @param draws     The number of stochastic parameter draws you would like to pull   
  #' output: list with site name, urban/rural grouping, iso code, and additional 'stoch_draws'
  #'         list with stochastic draws equal to the number of draws requested.
  
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


# transforming inputs ----------------------------------------------------------


time_transform <- function(x, time_divisor = 1, baseline_t = 0){
  #' Create a new time column for aggregation
  #'
  #' @param x Input data.frame
  #' @param time_divisor Aggregation level. Default = 1 (no aggregation).
  #' Setting to 365 would allow annual aggregation, 30 monthly, 7 weekly etc.
  #' @param baseline_t A baseline time to add to the time output
  #'
  #' @export
  #' 
  if(max(x$timestep) %% time_divisor != 0){
    warning("Number of timesteps not divisible exactly by level, group may be unequal")
  }
  
  x <- x |>
    dplyr::mutate(t = as.integer(ceiling(.data$timestep / time_divisor) + baseline_t)) |>
    dplyr::select(-"timestep")
  
  return(x)
}



drop_burnin <- function(x, burnin){
  #' Remove burn in period from simulation output
  #'
  #' @param x Input data.frame
  #' @param burnin Length of burn in period in days
  #'
  #' @export
  if(burnin >= max(x$timestep)){
    stop("burn in period must be < the maximum timestep")
  }
  if(burnin < 0){
    stop("burn in period must be positive")
  }
  
  
  x <- x[x$timestep > burnin, ]
  x$timestep <- as.integer(x$timestep - burnin)
  return(x)
}


aggregate_outputs<- function(dt, interval, folder){
  
  #' aggregate incident cases based on a pre-determined time interval (expressed in days). 
  #' aggregated to the country level or site level.
  #' 
  #' @param dt             raw model output from malariasimulation package
  #' @param interval       time period you would like to calculate incidence over, expressed in days.
  #' @param folder         folder to save aggregated output to
  
  #' output: data table with summed cases and rates over specified time interval and location grouping.
  
  # reformat case outputs to long
  # need clinical cases, severe cases, population, number treated, and number detected
  
  site<- unique(dt$site_name)
  
  message(paste0('aggregating outputs for site', site))
  
  require(data.table)

  dt <- dt |> 
    select(t, 
           ft,
           iso, 
           site_name,
           urban_rural,
           run,
           contains("n_inc_clin"),   # clinical incidence
           contains("n_inc_severe"), # severe incidence
           contains("n_age")         # population
           )
  
  dt<- dt |>
    tidyr::pivot_longer(
      cols = -c("t", "ft", "iso", "site_name", "urban_rural", "run")
    ) |>
    dplyr::mutate(
      name = stringr:: str_remove(.data$name, "_inc")
    ) 
  
  
  # classify parameters for aggregation
  dt<- data.table(dt)
  dt[name %like% "clinical", group:= 'clinical_incidence']
  dt[name %like% "severe", group:= 'severe_incidence']
  dt[name %like% "age", group:= 'population']
  

  dt<- dt |>
    tidyr::separate(
      col = "name", 
      into = c(NA, NA, "age_lower", "age_upper"),
      sep = "_",
      convert = TRUE
    ) 


# create aggregates
# for clinical and severe incidence, sum incident cases over time period -------
# for population, round population over the time period
  dt<- data.table(dt)
  
  grouping<- c('t', 'iso', 'site_name', 'urban_rural', 'run', 'age_lower', 'age_upper', 'group')
  dt[group== 'clinical_incidence', val:= sum(value), by= grouping]
  dt[group== 'severe_incidence', val:= sum(value), by= grouping]
  dt[group== 'population', val:= round(mean(value)), by= grouping]
  
dt<- unique(dt, by= grouping)
  
 dt<-  dt|>
    tidyr::pivot_wider(
      id_cols = c("t", "ft", "iso", "age_lower", "age_upper", "site_name", "urban_rural", "run"),
      names_from = 'group',
      values_from = 'val'
    )

 #  save output to folder
  write_rds(dt, file= paste0(folder, 'aggregated_output_', site, '.RDS'))
  
  message('completed aggregation')
  
}

aggregate_further<- function(dt){
  #' aggregate outputs further to the country or site level
  #' aggregated to the country level or site level.
  #' 
  #' @param location if 'iso', aggregate based on ISO code. If 'site_name', aggregate to site level.
  
  
  dt<- data.table(dt)
  
  grouping<- c('t', 'iso', location, 'run', 'age_lower', 'age_upper', 'group')
  
  dt[clinical_incidence:= sum(clinical_incidence), by= grouping]
  dt[severe_incidence:= sum(value), by= grouping]
  dt[population:= sum(population), by= grouping]
  
  dt<- unique(dt, by= grouping)
  
  message('done aggregating')
  
}
# calculating outputs ----------------------------------------------------------


calculate_deaths_ylls<- function(dt, cfr= 0.215, treatment_scaler= 0.5, lifespan= 63){
  
  #' Calculate deaths + years of life lost (YLLs) per GTS method.
  #' Where severe cases= 0, deaths= 0. Additionally remove a proportion of the cases that have received treatment.
  #' 
  #' @param dt               malariasimulation model outputs with columns 'severe' for severe incidence and 'n_treated' for number of individuals treated
  #' @param cfr              per GTS method, deaths are calculated using a case fatality ratio (CFR) value that is applied to severe incidence.
  #'                         see World Malaria Report and Wilson et al. for more information.
  #' @param scaler we remove a proportion of cases that have received treatment, assuming that 50% of treated cases remit and
  #'                         are no longer susceptible to mortality.
  #' @param lifespan         expected lifespan used to calculate YLLs. 
  #'                         YLLs are calculated by multiplying deaths by the number of years an individual was expected to live past their year of death.
  #'                         when comparing YLLs across different locations, it is recommended to use the same lifespan across YLL calculations.
  #'                         Typically, you should use the highest observed life expectancy in the region/ location you are studying.
  #' Output: data table with columns titled 'deaths' and 'yll'
  
  
  message('calculating deaths and YLLs')
  
  # add mortality rate
  dt<- dt |>
    mutate(mortality_rate = scaler * .data$severe_incidence)
  
  dt<- dt |>
    mutate(deaths =  .data$mortality_rate - (ft * treatment_scaler))
  
  
  dt<- data.table(dt)
  
  dt[deaths< 0 , deaths:= 0]
  dt[, deaths:= as.integer(deaths)]
  
  # calculate ylls
  # transform age to years
  dt[, age_years_lower:= age_lower / 365]
  dt[, age_years_upper:= age_upper / 365]
  
  dt <- dt |>
    mutate(yll= deaths * (lifespan - (age_years_upper - age_years_lower)/2))  
  
  message('completed calculation of deaths and YLLs')
  
  return(dt)
}

calculate_ylds_dalys<- function(dt, 
                                mild_dw= 0.006, 
                                moderate_dw= 0.051, 
                                severe_dw= 0.133,
                                clin_episode_length= 0.01375,
                                severe_episode_length= 0.04795,
                                interval= 365){
  
  #' Calculate Years Lived with Disability (YLDs) and Disability-Adjusted Life-Years
  #' based on disability weights from the Global Burden of Disease study.
  #' 
  #' Disability weights sourced here:
  #' https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fghdx.healthdata.org%2Fsites%2Fdefault%2Ffiles%2Frecord-attached-files%2FIHME_GBD_2017_DISABILITY_WEIGHTS_Y2018M11D08.XLSX&wdOrigin=BROWSELINK
  #' 
  #' Keep in mind this is an approximation of YLD estimation from the GBD study; disability due to comorbid conditions
  #' such as motor impairment and anemia are excluded.
  #' 
  #' For now, we assume that cases under 5 are moderate, due to the higher severity of malaria at younger ages.
  #' This assumption may be revisited in the future, 
  #' @param dt                    input dataset with columns 'clinical' for clinical incidence,
  #'                              'severe' for severe incidence, 'deaths', 'ylls', 'age_years_start', and 'age_years_end'
  #' @param mild_dw               disability weight for mild malaria
  #' @param moderate_dw           disability weight for moderate malaria
  #' @param severe_dw             disability weight for severe malaria
  #' @param clin_episode_length   length of an episode of clinical malaria
  #' @param severe_episode_length length of an episode of severe malaria 
  #' @param interval              period of time you are calculating YLDs for. This determines disease episode length.
  #' 
  #' Output: data table with columns 'ylds' and 'dalys'
  
  require(data.table)
  
  message('calculating YLDs and DALYs')
  
  # calculate YLDs  ---
  dt[ age_years_end < 5, 
      yld:= severe_incidence * severe_dw * severe_episode_length * interval + 
        clinical_incidence * moderate_dw * clin_episode_length * interval]
  
  dt[ age_years_end >= 5, 
      yld:= severe_incidence * severe_dw * severe_episode_length * interval + 
        clinical_incidence * mild_dw * clin_episode_length * interval]
  
  # calculate DALYs ---
  dt<- dt |>
    mutate(daly= yll+ yld)
  
  message('completed calculation of YLDs and DALYs')
  
  return(dt)
}


reformat_vimc_outputs<- function(dt){
  #' Reformat outputs for submission to VIMC
  #' @param dt data table you would like to reformat
  
  dt<- dt |>
    mutate(disease = 'Malaria',
           country_name = iso,
           country = iso) |>
    rename(cohort_size = population,
           cases = clinical_incidence,
           dalys = daly,
           age = age_years_start,
           year = t) |>
    select(disease, year, age, country, country_name, cohort_size, cases, dalys, deaths)
  
  return(dt)
}


calculate_cases_deaths_averted<- function(intvn_dt, bl_dt){
  #' Calculates cases, deaths, and dalys averted between intervention and baseline scenario.
  #' @param intvn_dt               model outputs for intervention scenario
  #' @param baseline               model outputs for baseline scenario
  #' 
  
  intvn_dt<- intvn_dt |>
    select(time, iso, run, age_years_start, age_years_end,
           clinical, severe, deaths, daly) |>
    rename(clinical_intvn= clinical,
           severe_intvn= severe,
           deaths_intvn= deaths,
           daly_intvn= daly)
  
  bl_dt<- bl_dt |>
    select(time, iso, run, age_years_start, age_years_end,
           clinical, severe, deaths, daly) |>
    rename(clinical_bl= clinical,
           severe_bl= severe,
           deaths_bl= deaths,
           daly_bl= daly)
  
  output<- left_join(bl_dt, intvn_dt, 
                     by= c('time', 'iso', 'age_years_start', 'age_years_end'))
  
  # calculate averted cases/deaths/dalys
  
  output<- output|>
    mutate(cases_averted = clinical_bl - clinical_intvn,
           severe_cases_averted = severe_bl - severe_intvn,
           deaths_averted = deaths_bl - deaths_intvn,
           dalys_averted = daly_bl - daly_intvn)
  
  message('calculated cases deaths dalys averted')
  
  return(output)
  
}






#' Extract basic rates from model output
#'
#' @param x Input data.frame
#' @inheritParams prevalence_format
#' @inheritParams time_transform
#'
#' @export
get_prevalence <- function(x, time_divisor, baseline_t, age_divisor){
  prevalence <- x |>
    prevalence_estimate() |>
    time_transform(time_divisor = time_divisor, baseline_t = baseline_t) |>
    prevalence_time_aggregate() |>
    prevalence_format(age_divisor = age_divisor)
  return(prevalence)
}




