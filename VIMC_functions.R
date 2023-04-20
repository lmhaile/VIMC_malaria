######################################################
## title    VIMC functions
## author   Lydia Haile
## purpose  functions used for VIMC malaria model runs
######################################################

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
      #min_ages = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95) *365,
      eir= site$eir$eir[1],
      overrides = list(human_population= 10000) # what size population is appropriate?
    )
    
    year<- 365
    # Set clinical incidence rendering
    params$clinical_incidence_rendering_min_ages = c(seq(0, 100, by = 1))*year
    params$clinical_incidence_rendering_max_ages = c(seq(0, 100, by = 1))*year
    
    # Set severe incidence rendering
    params$severe_incidence_rendering_min_ages = c(seq(1, 100, by = 1))*year
    params$severe_incidence_rendering_max_ages = c(seq(1, 100, by = 1))*year
    
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


aggregate_outputs<- function(dt, interval, sum_to_country){
  
  #' aggregate incident cases based on a pre-determined time interval (expressed in days). 
  #' aggregated to the country level or site level.
  #' 
  #' @param dt             raw model output from malariasimulation package
  #' @param interval       time period you would like to calculate incidence over, expressed in days.
  #' @param sum_to_country set to TRUE if you would like to sum values up to the country level (based on ISO code).
  #'                       if set to FALSE, outputs are aggregated to site level.
  #' output: data table with summed cases and rates over specified time interval and location grouping.
  
  # reformat case outputs to long
  # need clinical cases, severe cases, population, number treated, and number detected
  
  message(paste0('aggregating outputs by time interval: ', interval, ' days'))
  

  dt <- dt |> 
    select(timestep, 
           iso, 
           site_name,
           urban_rural,
           run,
           contains("n_inc_clin"), 
           contains("n_inc_sev"), 
           contains("n_age"), 
           contains('n_treated'), 
           contains('n_detect'))
  

  dt<- dt |>
    pivot_longer(c(contains("n_inc_clin"), contains("n_inc_sev"), contains("n_age")),
                 names_to = "age", 
                 values_to = "value") |>
    mutate(type = ifelse(grepl('n_inc_clinical', age), 'clinical', 'severe')) |>
    mutate(type = ifelse(grepl('n_age', age), 'population', type))|>
    mutate(age = gsub('n_inc_clinical_', '', age),
           age = gsub('n_inc_severe_', '', age),
           age = gsub('n_age_', '', age)) |>
    spread(key= type, value= value) |>
    separate(age, into = c("age_days_start", "age_days_end"), sep = "_") 
  
  
  # transform time interval
  dt <- dt |>
    mutate(time = as.integer(timestep/ interval))
  
  
  # aggregate outputs based on this interval, also by country
  
  if(sum_to_country== TRUE){
    
    dt<- dt |> 
      group_by(age_days_start, time, iso) |> 
      mutate(clinical = sum(clinical),
             severe = sum(severe),
             n_treated = sum(n_treated))
    
    # population is an average of summed counts over the time interval of interest
    
    dt<- dt |>
      group_by(age_days_start, timestep, iso) |>
      mutate(population = sum(population)) |>
      group_by(age_days_start, time, iso) |>
      mutate(population= round(mean(population)))
    
  } else {
    
    dt<- dt |> 
      group_by(age_days_start, time, site_name) |> 
      mutate(clinical = sum(clinical),
             severe = sum(severe),
             n_treated = sum(n_treated))
    
    # population is an average of summed counts over the time interval of interest
    
    dt<- dt |>
      group_by(age_days_start, timestep, site_name) |>
      mutate(population = sum(population)) |>
      group_by(age_days_start, time, site_name) |>
      mutate(population= round(mean(population)))
    
  }
    
    if(sum_to_country== TRUE){
      
    dt<- dt |>
      select(-timestep, -n_detect_365_36499, -n_detect_730_3649, -site_name, -urban_rural) |>
        distinct()
      
    } else {
      
      dt<- dt |>
        select(-timestep, -n_detect_365_36499, -n_detect_730_3649) |>
        distinct()
      
      
    }
  
  
  
  # calculate rates based on this interval
  dt<- dt |> 
    mutate(clin_rate = clinical/ population,
           severe_rate = severe/ population)

  return(dt)
  message('completed aggregation')
  
}


calculate_deaths_ylls<- function(dt, cfr= 0.215, treatment_scaler= 0.5, lifespan= 63){
  
  #' Calculate deaths + years of life lost (YLLs) per GTS method.
  #' Where severe cases= 0, deaths= 0. Additionally remove a proportion of the cases that have received treatment.
  #' 
  #' @param dt               malariasimulation model outputs with columns 'severe' for severe incidence and 'n_treated' for number of individuals treated
  #' @param cfr              per GTS method, deaths are calculated using a case fatality ratio (CFR) value that is applied to severe incidence.
  #'                         see World Malaria Report and Wilson et al. for more information.
  #' @param treatment_scaler we remove a proportion of cases that have received treatment, assuming that 50% of treated cases remit and
  #'                         are no longer susceptible to mortality.
  #' @param lifespan         expected lifespan used to calculate YLLs. 
  #'                         YLLs are calculated by multiplying deaths by the number of years an individual was expected to live past their year of death.
  #'                         when comparing YLLs across different locations, it is recommended to use the same lifespan across YLL calculations.
  #'                         Typically, you should use the highest observed life expectancy in the region/ location you are studying.
  #' Output: data table with columns titled 'deaths' and 'yll'
  
  
  message('calculating deaths and YLLs')
  
  # transform age into years
  dt<- dt|> 
    mutate(age_days_start= as.numeric(age_days_start),
           age_days_end= as.numeric(age_days_end))
  
  dt<- dt |>
    mutate(age_years_start= as.integer(age_days_start/ 365),
           age_years_end= as.integer(age_days_end/ 365))
  dt<- dt |>
    mutate(deaths =  cfr * severe - (n_treated * treatment_scaler))
  
  
  dt<- data.table(dt)
  
  dt[deaths< 0 , deaths:= 0]
  dt[, deaths:= as.integer(deaths)]
  
  # calculate ylls
  
  dt <- dt |>
    mutate(yll= deaths * (lifespan - (age_years_end - age_years_start)/2))  
  
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
      yld:= severe * severe_dw * severe_episode_length * interval + 
        clinical * moderate_dw * clin_episode_length * interval]
  
  dt[ age_years_end >= 5, 
      yld:= severe * severe_dw * severe_episode_length * interval + 
        clinical * mild_dw * clin_episode_length * interval]
  
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
           cases = clinical,
           dalys = daly,
           age = age_years_start,
           year = time) |>
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

