#############################################################
##  title   01_model_launch
##  author  Lydia Haile
##  purpose parameterize models and launch on HPC cluster
#############################################################

rm(list= ls())

# packages  --------------------------------------------------------------------
#install.packages('Q:/site_0.2.2.tar.gz')
library(scene)
library(tidyverse)
library(furrr)
library(data.table)
library(drat)
library(foresite)
library(dplyr)
library(mlgts)
library(tidyverse)
library(furrr)
library(scene)
library(site)
library(openxlsx)
library(wesanderson)
library(extrafont)
library(malariasimulation)
drat::addRepo("malariaverse", "file:\\\\projects.dide.ic.ac.uk/malaria/malariaverse/drat")


# custom functions to source  --------------------------------------------------
source("Q:/VIMC_malaria/functions/parameterizing_functions.R")
source("Q:/VIMC_malaria/functions/postprocessing_functions.R")
source("Q:/VIMC_malaria/functions/modelling_functions.R")


# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/'      #  directory where code is stored
malaria_dir<- 'Q:/VIMC/central_estimates/'      #  project directory where files are stored
input_dir<- 'Q:/VIMC'
setwd('Q:/')

# load in inputs from VIMC  ----------------------------------------------------
pop<- read.csv(paste0(input_dir, '/inputs/202212rfp-1_dds-202208_int_pop_both.csv'))
total_pop<- read.csv(paste0(input_dir, '/inputs/202212rfp-1_dds-202208_tot_pop_both.csv'))
le<- read.csv(paste0(input_dir, '/inputs/202212rfp-1_dds-202208_life_ex_both.csv'))
coverage<- read.csv(paste0(input_dir, '/inputs/coverage_202212rfp-1_malaria-mal4-default.csv'))
mort<- read.csv(paste0(input_dir, '/inputs/mortality.csv'))

# format mortality data --------------------------------------------------------
mort<- mort[, c('age_to', 'year', 'value')]
mort<- data.table(mort)
mort[age_to== 0, age_to:= 1]
mort[, age_to:= age_to * 365]

# cut off mortality data before 2000
mort<- mort[year >= 2000 & year <= 2050]

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

# pull site data  --------------------------------------------------------------
site_data<- foresite::NGA
#site_data<- site::single_site(site_data, 50)

# plot initial vaccine coverage  -----------------------------------------------
plot_interventions_combined(
  interventions = site_data$interventions,
  population = site_data$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment","SMC", "PMC")
)

#' Format input parameters for VIMC sites
#'
#' @param   site_data          site data file for a particular country
#' @param   demography         if TRUE, incorporate VIMC demography inputs
#' @param   population         human population to run model on
#' @param   scenario           'baseline' or 'intervention'
#' @param   short_run          if true, shorten site file data to year pre-2020
#' @param   min_ages           minimum ages (in timesteps) for incidence/ prevalence outputs
#' @param   max_ages           maximum ages (in timesteps) for incidence/ prevalence outputs
#' @param   tag                description of model runs you would like to carry out
#' @returns list with input parameters and demographic info, as inputs for run_malaria_model


prep_model_launch<- function(site_data, 
                             vimc_mortality,
                             vimc_population,
                             population, 
                             scenario, 
                             min_ages, 
                             max_ages, 
                             tag){
  
  isos<- unique(site_data$sites$iso3c)
  
  prep_site_data<- function(num){
  
  message(paste0('prepping site ', num))
  site<- site::single_site(site_file= site_data, index= num) 
  
  # get site info
  site_name<- site$sites$name_1
  ur<- site$sites$urban_rural
  iso<- site$sites$iso3c
  
  # create a directory to save your output
  if(dir.exists(paste0('Q:/VIMC/central_estimates/input_parameters/', iso))== FALSE){
    dir.create(paste0('Q:/VIMC/central_estimates/input_parameters/', iso))
  }
  
  if(dir.exists(paste0('Q:/VIMC/central_estimates/input_parameters/', iso, '/', tag))== FALSE){
    dir.create(paste0('Q:/VIMC/central_estimates/input_parameters/', iso, '/', tag))
  }

  if(vimc_mortality== TRUE){
  # bind on a year for youngest age group
  youngest<- data.table(site$demography)[age_upper== min(site$demography$age_upper)]
  mort_dt<- rbind(youngest, mort)
  site$demography<- mort_dt
  }
  
  if (vimc_population== TRUE){
  # replace population in site file with population from VIMC inputs  ------------
  site$population<- merge(site$population, total_pop[, c('year', 'value')], by = 'year')
  site$population<- site$population |>
    select(-pop)
  site$population<- site$population |>
    rename(pop = value)
  }
  message(paste0('prepping inputs for site ', site_name, ' ', ur))
  
  if (scenario== 'baseline'){
    
    # # expand scenario out to 2050
    site <- set_vaccine_coverage(
      site,
      change = FALSE,
      terminal_year = 2050,
      rtss_target = NULL,
      rtss_year = NULL
    )

  }else if (scenario== 'intervention') {
    site <- set_vaccine_coverage(
      site,
      change = TRUE,
      terminal_year = 2050,
      rtss_target = 0.8,
      rtss_year = 2023
    )
  }
  
 # calibration<- quick_calibration(site)
  
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
  
  inputs<- list('param_list'= params, 
                'site_name'= site_name, 
                'ur'= ur, 
                'iso'= iso,
                'scenario'= scenario,
                'tag'= tag)
  
  write_rds(
    inputs,
    paste0('Q:/VIMC/central_estimates/input_parameters/',iso, '/', tag, '/', site_name, '_', ur, '_', scenario,
      '.rds'
    )
  )
  
  }
  lapply(seq(1:nrow(site_data$sites)), prep_site_data)
  return(message(paste0('prepped outputs for model run: ', tag, ', country: ', isos)))
  
}

year<- 365
model_input<- prep_model_launch(site_data, 
                                scenario= 'baseline', 
                                vimc_mortality= TRUE,
                                vimc_population= FALSE,
                                min_ages = c(seq(0, 14, by= 1), seq(15, 80, by= 15)) * year,
                                max_ages = c(seq(1, 15, by= 1), seq(30, 95, by= 15)) * year -1,
                                population = 50000,
                                tag= 'mortality_change_test_run')






# submit jobs to cluster  ------------------------------------------------------
# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation")

# save a context ---------------------------------------------------------------
# additional script contains helper functions for larger scale model runs
ctx <- context::context_save(
  'pkgs',
  packages = packages,
  package_sources = src,
  sources = 'Q:/VIMC_malaria/functions/modelling_functions.R'
)

config <- didehpc::didehpc_config(
  use_rrq = TRUE,
  cores = 1,
  cluster = "fi--didemrchnb" , 
  template = "GeneralNodes",
  shares = didehpc::path_mapping('nas2', #
                                 path_local = 'M:/', 
                                 path_remote = '\\\\fi--didenas1.dide.ic.ac.uk', 
                                 drive_remote = 'M:')
)

# load context into queue
obj <- didehpc::queue_didehpc(ctx, config)
#obj <- didehpc::queue_didehpc(ctx)

# launch jobs ---------------------------------------------
iso<- 'NGA'
tag<- 'mortality_change_test_run'

filepaths<- list.files(paste0(malaria_dir, 'input_parameters/', iso, '/', tag),
                       full.names= T)

jobs<- obj$lapply(filepaths, run_malaria_model)



# test if the jobs begin before launching on cluster
filepath<- filepaths[[1]]
lapply(filepath, run_malaria_model)
