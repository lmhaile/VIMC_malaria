#############################################################
##  title   01_model_launch
##  author  Lydia Haile
##  purpose parameterize models and launch on HPC cluster
#############################################################

rm(list= ls())

# packages  --------------------------------------------------------------------
remotes::install_github('mrc-ide/site')
install.packages('Q:/site_0.2.2.tar.gz')
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

# merge on mortality inputs ----------------------------------------------------
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

# replace demography input with VIMC input
site$demography<- mort[year < 2051]

# replace population in site file with population from VIMC inputs  ------------
site$population<- merge(site$population, total_pop[, c('year', 'value')], by = 'year')
site$population<- site$population |>
select(-pop)
site$population<- site$population |>
rename(pop = value)

# plot initial vaccine coverage  -----------------------------------------------
plot_interventions_combined(
  interventions = site$interventions,
  population = site$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment","SMC", "PMC")
)

# make copies of site data (one for baseline scenario and one for intervention)
baseline<- copy(site)
intvn <- copy(site)

# set vaccine coverage ---------------------------------------------------------
# run a basic site to see if cases + deaths match what is expected
# intvn <- set_vaccine_coverage(
#   intvn,
#   change = TRUE,
#   terminal_year = 2050,
#   rtss_target = 0.8,
#   rtss_year = 2023
# )
# 
# baseline <- set_vaccine_coverage(
#   baseline,
#   change = FALSE,
#   terminal_year = 2050,
#   rtss_target = 0.8,
#   rtss_year = 2023
# )

# plot the changes you made ----------------------------------------------------
plot_interventions_combined(
  interventions = intvn$interventions,
  population = intvn$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment","SMC", "PMC")
)


# plot baseline to make sure they look different  ------------------------------
plot_interventions_combined(
  interventions = baseline$interventions,
  population = baseline$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment","SMC", "PMC")
)

year<- 365
# prep central burden estimate inputs ------------------------------------------
bl <- prep_inputs(
  baseline,
  folder = paste0(malaria_dir, '/central_estimates/baseline/'),
  population = 100,
  min_ages = 0 * year,
  max_ages = 15 * year -1
)

# int<- prep_inputs(intvn, 
#                   folder= paste0(malaria_dir, '/central_estimates/intervention/'),
#                   population = 100,
#                   min_ages= c(seq(0, 80, by = 20)) * year,
#                   max_ages = c(seq(20, 100, by = 20)) * year -1)


# prep stochastic burden estimate inputs
#int_stochastic<- lapply(int, prep_stochastic_inputs, draws= 10)
#int_stochastic<- flatten(int_stochastic)

#bl_stochastic<- lapply(bl, prep_stochastic_inputs, draws= 10)
#bl_stochastic<- flatten(bl_stochastic)

# submit jobs to cluster  ------------------------------------------------------
message(paste0('submitting ', length(bl),  ' central burden jobs'))
#message(paste0('submitting ', length(int),  ' central burden jobs'))


# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation")

# save a context (working environment for your code) ---------------------------
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
  cluster = "fi--didemrchnb" ,
  #"fi--dideclusthn", # , "fi--didemrchnb""fi--didemrchnb"
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

# run central burden baseline jobs ---------------------------------------------
central_fold<- paste0(malaria_dir, 
                      '/central_estimates/baseline/') # folder you would like to save outputs in
stochastic_fold<-  paste0(malaria_dir, 
                          '/stochastic_estimates/baseline/')

central_baseline_jobs <- obj$lapply(bl, 
                                    run_malaria_model, 
                                    folder= central_fold,
                                    stochastic_run= FALSE,
                                    tag= 'no_changes')


#stochastic_baseline_jobs<- obj$lapply(bl_stochastic, 
#run_malaria_model, 
#folder= stochastic_fold,
#stochastic_run= T)


# run central burden  intervention jobs ----------------------------------------
central_fold<- paste0(malaria_dir, 
                      '/central_estimates/intervention/') # folder you would like to save outputs in
stochastic_fold<-  paste0(malaria_dir, 
                          '/stochastic_estimates/intervention/')

# central_intvn_jobs <- obj$lapply(int, 
#                                  run_malaria_model, 
#                                  folder= central_fold,
#                                  stochastic_run= FALSE,
#                                  tag= 'VIMC_mortality')

#stochastic_baseline_jobs<- obj$lapply(int_stochastic,  
#run_malaria_model, 
#folder= stochastic_fold)
test<- lapply(int,
              run_malaria_model, 
              folder= central_fold)

