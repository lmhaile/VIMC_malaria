#   archive functions ----------------------------------------------------------

#' aggregate incident cases based on a pre-determined time interval (expressed in days). 
#' aggregated to the country level or site level.
#' 
#' @param dt             raw model output from malariasimulation package
#' @param interval       time period you would like to calculate incidence over, expressed in days.
#' @param folder         folder to save aggregated output to

#' output: data table with summed cases and rates over specified time interval and location grouping.

# reformat case outputs to long
# need clinical cases, severe cases, population, number treated, and number detected



aggregate_outputs<- function(dt, interval, folder){
  
  dt<- readRDS(filepath)
  
  site<- unique(dt$site_name)
  #dt <- time_transform(x= dt, time_divisor = 365, baseline_t = 0) # transform time into annual outputs
  
  message(paste0('aggregating outputs for site ', site))
  
  require(data.table)
  
  dt <- dt |> 
    select(timestep, 
           ft,
           iso, 
           site_name,
           urban_rural,
           scenario,
           description,
           contains("n_inc_clin"),   # clinical incidence
           contains("n_inc_severe"), # severe incidence
           contains("n_age")         # population
    )
  
  dt<- dt |>
    tidyr::pivot_longer(
      cols = -c("timestep", "ft", "iso", "site_name", "urban_rural", "scenario", 'description')
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
  
  dt<- data.table(dt)
  
  # transform age to years
  dt[, age_years_lower:= round(age_lower / 365)]
  dt[, age_years_upper:= round(age_upper / 365)]
  
  
  # create aggregates
  # for clinical and severe incidence, sum incident cases over time period -------
  # for population, round population over the time period
  dt<- data.table(dt)
  
  # grouping<- c('t', 'iso', 'site_name', 'urban_rural', 'scenario', 'description', 'age_lower', 'age_upper', 'group')
  # dt[group== 'clinical_incidence', val:= sum(value), by= grouping]
  # dt[group== 'severe_incidence', val:= sum(value), by= grouping]
  # dt[group== 'population', val:= round(mean(value)), by= grouping]
  # 
  dt<- unique(dt, by= grouping)
  
  dt<-  dt|>
    tidyr::pivot_wider(
      id_cols = c("timestep", "ft", "iso", "age_years_lower", "age_years_upper", "site_name", "urban_rural", "scenario", "description"),
      names_from = 'group',
      values_from = 'value'
    )
  
  return(dt)
  
  message('completed aggregation')
  
}



# calculating outputs ----------------------------------------------------------

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

calculate_deaths_ylls<- function(dt, cfr= 0.215, scaler= 0.57, lifespan= 63){
  
  message('calculating deaths and YLLs')
  # add mortality rate
  dt<- dt |>
    mutate(mortality_rate = scaler * .data$severe_incidence)
  
  dt<- dt |>
    mutate(deaths =  .data$mortality_rate - (ft * scaler))
  
  
  dt<- data.table(dt)
  
  dt[deaths< 0 , deaths:= 0]
  dt[, deaths:= as.integer(deaths)]
  
  # calculate ylls
  dt <- dt |>
    mutate(yll= round(deaths * (lifespan - (age_years_upper - age_years_lower)/2)))
  
  message('completed calculation of deaths and YLLs')
  
  return(dt)
}


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
#' @returns data table with columns 'ylds' and 'dalys'

calculate_ylds_dalys<- function(dt, 
                                mild_dw= 0.006, 
                                moderate_dw= 0.051, 
                                severe_dw= 0.133,
                                clin_episode_length= 0.01375,
                                severe_episode_length= 0.04795,
                                interval= 365){
  
  
  require(data.table)
  
  message('calculating YLDs and DALYs')
  
  # calculate YLDs  ---
  dt[ age_years_upper < 5, 
      yld:= severe_incidence * severe_dw * severe_episode_length * 365 * interval + 
        clinical_incidence * moderate_dw * clin_episode_length * 365 * interval]
  
  dt[ age_years_upper >= 5, 
      yld:= severe_incidence * severe_dw * severe_episode_length * 365 * interval + 
        clinical_incidence * mild_dw * clin_episode_length * 365 * interval]
  
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
           country = iso,
           t = t + 1999) |>
    rename(cohort_size = population,
           cases = clinical_incidence,
           dalys = daly,
           age = age_years_lower,
           year = t,
           description = description) |>
    select(disease, year, age, country, country_name, site_name, urban_rural, scenario, description, cohort_size, cases, dalys, deaths)
  
  return(dt)
}


#' Calculates cases, deaths, and dalys averted between intervention and baseline scenario.
#' @param intvn_dt               model outputs for intervention scenario
#' @param baseline               model outputs for baseline scenario
#' 


calculate_cases_deaths_averted<- function(intvn_dt, bl_dt){
  
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










