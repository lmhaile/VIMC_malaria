#############################################################
##  title   01_model_launch
##  author  Lydia Haile
##  purpose parameterize models and launch on HPC cluster
#############################################################

rm(list= ls())

# packages  --------------------------------------------------------------------
#remotes::install_github('mrc-ide/site')
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
# install.packages("foresite", type = "source") # v0.1.0
# install.packages('site')
#install.packages('malariasimulation')
source("Q:/VIMC_malaria/VIMC_functions.R", echo=TRUE)
source("Q:/VIMC_malaria/run_malaria_model.R")

#remotes::install_github('mrc-ide/malariasimulation@dev')

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/' #  directory where code is stored
malaria_dir<- 'Q:/VIMC_files'      #  project directory where files are stored
setwd('Q:/')

# load in inputs from VIMC  ----------------------------------------------------
pop<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_int_pop_both.csv'))
total_pop<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_tot_pop_both.csv'))
le<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_life_ex_both.csv'))
coverage<- read.csv(paste0(malaria_dir, '/inputs/coverage_202212rfp-1_malaria-mal4-default.csv'))
mort<- read.csv(paste0(malaria_dir, '/inputs/mortality.csv'))
#fert<- read.csv(paste0(code_dir, '/inputs/fertility.csv'))

# pull site data  --------------------------------------------------------------
site<- foresite::NGA
site<- site::single_site(site, 50)

# plot initial vaccine coverage  -----------------------------------------------
plot_interventions_combined(
  interventions = site$interventions,
  population = site$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment","SMC", "PMC")
)


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
  site$population<- merge(site$population, total_pop[, c('year', 'value')], by = 'year')#

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
  sources = 'Q:/VIMC_malaria/run_malaria_model.R'
)

config <- didehpc::didehpc_config(
  use_rrq = FALSE,
  cores = 1,
  cluster = "fi--didemrchnb" , #"fi--dideclusthn", # , "fi--didemrchnb""fi--didemrchnb"
  template = "GeneralNodes"
  ## use for the wpia cluster
  #parallel = FALSE
)

# load context into queue
obj <- didehpc::queue_didehpc(ctx, config)

# obj <- didehpc::queue_didehpc(ctx, provision = "later")
# obj$install_packages("mrc-ide/individual")
# obj$install_packages("mrc-ide/malariaEquilibrium")
# obj$install_packages("mrc-ide/malariasimulation")

obj <- didehpc::queue_didehpc(ctx)

# run test job ---------------------------------------------
job<- obj$enqueue(run_malaria_model(model_input))

