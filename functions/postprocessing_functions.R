######################################################
## title    Post processing functions
## author   Lydia Haile
## purpose  functions used to process model outputs
##          portions pulled from scene package
######################################################

#' General VIMC outputs from raw model output
#' @param filepath filepath where raw model output is stored
#'  
#' @returns processed outputs for sites in output directory
generate_vimc_output<- function(filepath){
  
  dt<- readRDS(filepath)
  
  site_name<- unique(dt$site_name)
  ur<- unique(dt$ur)
  scenario<- unique(dt$scenario)
  description<- unique(dt$description)
  country<- unique(dt$iso)
  
  message(paste0('generating output for site ', site_name, ' ', ur))
  
  dt <- postie::get_rates(
    dt,
    time_divisor = 365,
    baseline_t = 2000,
    age_divisor = 365,
    scaler = 0.215,
    treatment_scaler = 0.42,
  ) 
  
  # merge in population from site files (as we do not have VIMC inputs for this yet)
  pop<- foresite::get_site(country)
  pop<- pop$population |>
    filter(name_1== site_name & urban_rural == ur ) |>
    select(year, pop) |>
    rename(t = year)
  
  dt<- merge(dt, pop, by = 't')
  
  dt<- dt |>
    mutate(disease = 'Malaria',
           cases = round(clinical * pop * prop_n),
           deaths = round(mortality * pop * prop_n),
           dalys = round(dalys_pp * pop * prop_n),
           population = round(pop * prop_n),
           country = country,
           country_name = countrycode::countrycode(sourcevar= country, origin= 'iso3c', destination= 'country.name'),
           site_name = site_name,
           urban_rural = ur,
           scenario = scenario,
           description = description) |>
    rename(year = t,
           age = age_lower,
           cohort_size = population) |>
    select(disease, year, age, country, country_name, site_name, urban_rural, prop_n,
           scenario, description, cohort_size, cases, dalys, deaths, clinical, mortality, dalys_pp) |>
    mutate(cases = if_else(is.na(cases), 0, cases),
           deaths = if_else(is.na(deaths), 0, deaths),
           dalys = if_else(is.na(dalys), 0, dalys),
           mortality= if_else(is.na(mortality), 0, mortality),
           clinical= if_else(is.na(clinical), 0, clinical),
           dalys= if_else(is.na(dalys), 0, dalys))  
  
  return(dt)
}

#' Aggregate VIMC output up to the country level 
#' @param filepath filepath where raw model output is stored
#'  
#' @returns aggregated country-level outputs
aggregate_output<- function(dt){
  
  dt <- dt |>
    group_by(year, age, country) |>
    mutate(cases = sum(cases),
           deaths = sum(deaths),
           dalys = sum(dalys),
           cohort_size = sum(cohort_size)) |>
    ungroup() |>
    select(-site_name, -urban_rural)
  
  
  dt<- distinct(dt, year, age, country, .keep_all = TRUE)

  return(dt)
}


#' remove extraneous columns for submission
#' @param dt VIMC output
#'  
#' @returns final data frame for submission
format_for_submission<- function(dt){
  
  dt<- dt |>
    select(disease, year, age, country, country_name, cohort_size, cases, dalys, deaths)
    
  return(dt)
}

