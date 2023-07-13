######################################################
## title    Post processing functions
## author   Lydia Haile
## purpose  functions used to process model outputs
##          portions pulled from scene package
######################################################

#' Create a new time column for aggregation
#'
#' @param x Input data.frame
#' @param time_divisor Aggregation level. Default = 1 (no aggregation).
#' Setting to 365 would allow annual aggregation, 30 monthly, 7 weekly etc.
#' @param baseline_t A baseline time to add to the time output
#'
#' @export
#' 

time_transform <- function(x, time_divisor = 1, baseline_t = 0){
  if(max(x$timestep) %% time_divisor != 0){
    warning("Number of timesteps not divisible exactly by level, group may be unequal")
  }
  
  x <- x |>
    dplyr::mutate(t = as.integer(ceiling(.data$timestep / time_divisor) + baseline_t)) |>
    dplyr::select(-"timestep")
  
  return(x)
}


#' Remove burn in period from simulation output
#'
#' @param x Input data.frame
#' @param burnin Length of burn in period in days
#'
#' @export

drop_burnin <- function(x, burnin){
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
  
  
  site<- unique(dt$site_name)
  
  message(paste0('aggregating outputs for site ', site))
  
  require(data.table)
  
  dt <- dt |> 
    select(t, 
           ft,
           iso, 
           site_name,
           urban_rural,
           scenario,
           tag,
           contains("n_inc_clin"),   # clinical incidence
           contains("n_inc_severe"), # severe incidence
           contains("n_age")         # population
    )
  
  dt<- dt |>
    tidyr::pivot_longer(
      cols = -c("t", "ft", "iso", "site_name", "urban_rural", "scenario", 'tag')
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
  
  grouping<- c('t', 'iso', 'site_name', 'urban_rural', 'scenario', 'tag', 'age_lower', 'age_upper', 'group')
  dt[group== 'clinical_incidence', val:= sum(value), by= grouping]
  dt[group== 'severe_incidence', val:= sum(value), by= grouping]
  dt[group== 'population', val:= round(mean(value)), by= grouping]
  
  dt<- unique(dt, by= grouping)
  
  dt<-  dt|>
    tidyr::pivot_wider(
      id_cols = c("t", "ft", "iso", "age_years_lower", "age_years_upper", "site_name", "urban_rural", "scenario", "tag"),
      names_from = 'group',
      values_from = 'val'
    )
  
  return(dt)
  
  message('completed aggregation')
  
}


#' aggregate outputs further to the country or site level
#' aggregated to the country level or site level.
#' 
#' @param location if 'iso', aggregate based on ISO code. If 'site_name', aggregate to site level.


aggregate_further<- function(dt){
  
  dt<- data.table(dt)
  
  grouping<- c('t', 'iso', location, 'run', 'age_years_lower', 'age_years_upper', 'group')
  
  dt[clinical_incidence:= sum(clinical_incidence), by= grouping]
  dt[severe_incidence:= sum(value), by= grouping]
  dt[population:= sum(population), by= grouping]
  
  dt<- unique(dt, by= grouping)
  
  message('done aggregating')
  
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
           description = tag) |>
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





#' Extract basic rates from model output
#'
#' @param x Input data.frame
#' @inheritParams rates_format
#' @inheritParams time_transform
#' @inheritParams mortality_rate
#' @inheritParams treatment_scaling
#' @param aggregate_age Aggregate output over age groups
#'
#' @export
get_rates <- function(x, time_divisor, baseline_t, age_divisor, scaler, treatment_scaler, baseline_treatment, aggregate_age = FALSE){
  rates <- x |>
    rates_column_check() |>
    rates_transform() |>
    time_transform(time_divisor = time_divisor, baseline_t = baseline_t) |>
    rates_time_aggregate() |>
    rates_format(age_divisor = age_divisor) |>
    treatment_scaling(treatment_scaler = treatment_scaler, baseline_treatment = baseline_treatment) |>
    mortality_rate(scaler = scaler)
  
  if(aggregate_age){
    rates <- rates |>
      rates_age_aggregate()
  }
  return(rates)
}

#' Rates input checks
#'
#' Checks that the required columns are present and that the age-ranges for
#' required columns align.
#'
#' @param x Input data.frame
rates_column_check <- function(x){
  cols <- colnames(x)
  if(!"timestep" %in% cols){
    stop("required column `timestep` missing")
  }
  if(!"ft" %in% cols){
    stop("required column `ft` missing")
  }
  if(sum(grepl("n_inc_clinical", cols)) == 0){
    stop("required columns `n_inc_clinical_...` missing")
  }
  if(sum(grepl("n_inc_severe", cols)) == 0){
    stop("required columns `n_inc_severe_...` missing")
  }
  if(sum(grepl("n_age_", cols)) == 0){
    stop("required columns `n_age_...` missing")
  }
  clinical_cols <- cols[grepl("n_inc_clinical", cols)]
  clinical_age_ranges <- stringr::str_split(stringr::str_replace(clinical_cols, "n_inc_clinical_", ""), "_")
  severe_cols <- cols[grepl("n_inc_severe", cols)]
  severe_age_ranges <- stringr::str_split(stringr::str_replace(severe_cols, "n_inc_severe_", ""), "_")
  n_age_cols <- cols[grepl("n_age", cols)]
  n_age_ranges <- stringr::str_split(stringr::str_replace(n_age_cols, "n_age_", ""), "_")
  if(!identical(clinical_age_ranges, severe_age_ranges) | !identical(clinical_age_ranges, n_age_ranges)){
    stop("Age ranges for `n_inc_clinical_...` and `n_inc_severe_...` and `n_age_...` outputs must be the same")
  }
  return(x)
}

#' Transform rates into long form
#'
#' @param x Input data.frame
rates_transform <- function(x){
  x <- x |>
    dplyr::select(
      "timestep",
      "ft",
      dplyr::contains("n_age"),
      dplyr::contains("n_inc_clinical"),
      dplyr::contains("n_inc_severe"),
    ) |>
    tidyr::pivot_longer(
      cols = -c("timestep", "ft")
    ) |>
    dplyr::mutate(
      name = stringr::str_remove(.data$name, "_inc")
    ) |>
    tidyr::separate(
      col = "name",
      into = c(NA, "name", "age_lower", "age_upper"),
      sep = "_",
      convert = TRUE
    ) |>
    tidyr::pivot_wider(
      id_cols = c("timestep", "age_lower", "age_upper", "ft")
    )
  return(x)
}

#' Aggregate rates output over t
#'
#' @param x Input data.frame
rates_time_aggregate <- function(x){
  x <- x |>
    dplyr::summarise(
      dplyr::across(c("clinical", "severe"), sum),
      dplyr::across(c("age", "ft"), mean),
      .by = c("t", "age_lower", "age_upper")
    )
  return(x)
}


#' Format rates
#'
#' Create rates from counts and specifies age output units
#'
#' @param x Input data.frame
#' @param age_divisor Aggregation level. For example setting to 365 will return
#' age units in years
rates_format <- function(x, age_divisor = 365){
  if(age_divisor < 1){
    stop("age_divisor must be > 1")
  }
  
  x <- x |>
    dplyr::mutate(
      prop_age = .data$age / sum(.data$age),
      .by = "t"
    ) |>
    dplyr::mutate(
      clinical = .data$clinical / .data$age,
      severe = .data$severe / .data$age,
      age_lower = round(.data$age_lower / age_divisor),
      age_upper = round(.data$age_upper / age_divisor)) 
  #dplyr::select(-"age")
  return(x)
}

#' Aggregate rates output over age
#'
#' @param x Input data.frame
rates_age_aggregate <- function(x){
  x <- x |>
    dplyr::summarise(
      clinical = stats::weighted.mean(.data$clinical, .data$prop_age),
      severe  = stats::weighted.mean(.data$severe, .data$prop_age),
      mortality = stats::weighted.mean(.data$mortality, .data$prop_age),
      .by = "t"
    )
  return(x)
}

#' Adjust severe (and downstream mortality) rates as a result of treatment coverage.
#'
#' @param x Input data.frame
#' @param treatment_scaler The impact of treatment coverage on progression to severe disease and death.
#' A highly uncertain parameter, the analysis by
#'  \href{https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(15)00423-5/fulltext}{Griffin et al (2016)}
#' sampled from a uniform diastribution (0, 1).
#' @param baseline_treatment The proportion of uncomplicated malaria cases that are effectively treated historically.
treatment_scaling <- function(x, treatment_scaler, baseline_treatment = 0){
  if(treatment_scaler > 1 || treatment_scaler < 0){
    stop("treatment_scaler must be between 0 and 1")
  }
  if(baseline_treatment > 1 || baseline_treatment < 0){
    stop("treatment_scaler must be between 0 and 1")
  }
  ts <- (1 - treatment_scaler)
  x <- x |>
    dplyr::mutate(severe = .data$severe * ((ts * .data$ft + (1 - .data$ft)) / (ts * baseline_treatment + (1 - baseline_treatment)))) |>
    dplyr::select(-"ft")
  return(x)
}

mortality_rate <- function(x, scaler){
  if(scaler > 1 || scaler < 0){
    stop("scaler must be between 0 and 1")
  }
  
  x <- x |>
    dplyr::mutate(mortality = scaler * .data$severe)
  return(x)
}



