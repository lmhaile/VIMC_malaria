################################################################################
##  title   site file test case 
##  author  Lydia Haile
##  purpose test site file code locally to profile jobs
################################################################################

rm(list= ls())
# packages  --------------------------------------------------------------------
library(scene)
library(tidyverse)
library(furrr)
library(data.table)
library(drat)
library(foresite)
library(dplyr)
library(mlgts)
library(furrr)
library(site)
library(openxlsx)
library(wesanderson)
library(extrafont)
library(malariasimulation)

malaria_dir<- 'Q:/VIMC_files'      #  project directory where files are stored

# load in inputs from VIMC  ----------------------------------------------------
pop<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_int_pop_both.csv'))
total_pop<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_tot_pop_both.csv'))
le<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_life_ex_both.csv'))
coverage<- read.csv(paste0(malaria_dir, '/inputs/coverage_202212rfp-1_malaria-mal4-default.csv'))
mort<- read.csv(paste0(malaria_dir, '/inputs/mortality.csv'))
year<- 365
# objects ----------------------------------------------------------------------
min = c(seq(0, 14, by= 1), seq(15, 80, by= 15)) * year
max = c(seq(1, 15, by= 1), seq(35, 95, by= 15)) * year -1

# directory --------------------------------------------------------------------
folder<- 'Q:/VIMC_files/central_estimates/baseline/'

# pick one site
site<- foresite::NGA
site_data<- site::single_site(site, 50)
site<- site::single_site(site_file= site_data, index= 1) 

prep_model_launch<- function(site, population, scenario, min_ages= min, max_ages= max, tag= 'population_50k'){

mort<- read.csv(paste0(malaria_dir, '/inputs/mortality.csv'))

#get site info
site_name<- site$sites$name_1
ur<- site$sites$urban_rural
iso<- site$sites$iso3c


# format mortality data
mort<- mort[, c('age_to', 'year', 'value')]
mort<- data.table(mort)
mort[age_to== 0, age_to:= 1]
mort[, age_to:= age_to * 365]

# cut off mortality data before 2000
mort<- mort[year >= 2000]

# maybe attempt to bind mortality rate data directly into site file
# instead of using set_demography
mort<- mort |>
  rename(age_upper = age_to,
         mortality_rate = value) |>
  mutate(age_upper = age_upper/ 365) |>
  mutate(iso3c= 'NGA',
         country= 'Nigeria')

# align age groups
mort[age_upper > 1, age_upper := age_upper + 1]
mort[age_upper== 121, age_upper:= 200]

# reorder
mort<- mort |>
  select(iso3c, country, age_upper, year, mortality_rate)

# bind on a year for youngest age group
youngest<- data.table(site$demography)[age_upper== min(site$demography$age_upper)]
mort<- rbind(youngest, mort)

# replace population in site file with population from VIMC inputs  ------------
site$population<- merge(site$population, total_pop[, c('year', 'value')], by = 'year')
site$population<- site$population |>
select(-pop)
site$population<- site$population |>
rename(pop = value)

message(paste0('prepping inputs for site ', site_name, ' ', ur))

if (scenario== 'baseline'){
  
  # expand scenario out to 2050
  site <- set_vaccine_coverage(
    site,
    change = FALSE,
    terminal_year = 2050,
    rtss_target = 0.8,
    rtss_year = 2023
  ) 
  
}else if (scenario== 'intervention') {
  site <- set_vaccine_coverage(
    intvn,
    change = TRUE,
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
params$clinical_incidence_rendering_min_ages =  min_ages
params$clinical_incidence_rendering_max_ages = max_ages 

# Set age group rendering
params$age_group_rendering_min_ages = min_ages 
params$age_group_rendering_max_ages = max_ages 

inputs<- list('param_list'= params, 'site_name'= site_name, 'ur'= ur, 'iso'= iso, 'tag'= tag, 'scenario'= scenario)

return(inputs)

}


model_input<- prep_model_launch(site, scenario= 'baseline', population = 50000)

# run the model ----------------------------------------------------------------

run_malaria_model<- function(model_input) {
  
  params<-model_input$param_list
  params$progress_bar<- TRUE
  timesteps<<- params$timesteps
  
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

run_malaria_model(model_input)
