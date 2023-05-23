#############################################################
##  title   01_model_launch
##  author  Lydia Haile
##  purpose parameterize models and launch on HPC cluster
#############################################################

rm(list= ls())


# packages  --------------------------------------------------------------------
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
library(openxlsx)
library(wesanderson)
library(extrafont)
library(malariasimulation)
source("Q:/VIMC_malaria/VIMC_functions.R", echo=TRUE)

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/' #  directory where code is stored
malaria_dir<- 'Q:/VIMC_files'      #  project directory where files are stored
setwd('Q:/')

# load in inputs from VIMC  ----------------------------------------------------
pop<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_int_pop_both.csv'))
total_pop<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_tot_pop_both.csv'))

# le<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_life_ex_both.csv'))
coverage<- read.csv(paste0(malaria_dir, '/inputs/coverage_202212rfp-1_malaria-mal4-default.csv'))
mort<- read.csv(paste0(malaria_dir, '/inputs/mortality.csv'))
# fert<- read.csv(paste0(code_dir, '/inputs/fertility.csv'))

# pull site data  --------------------------------------------------------------
site<- foresite::MWI
site<- site::single_site(site, 54)


# format mortality data
mort<- mort[, c('age_to', 'year', 'value')]
mort<- data.table(mort)
mort[age_to== 0, age_to:= 1]
mort[, age_to:= age_to * 365]

# bind on a baseline year
bl_yr<- mort[year== min(mort$year)]
bl_yr[,year:= 0]
mort<- rbind(bl_yr, mort)

mort_mat<- dcast(mort, year ~ age_to, value.var = 'value')
setnames(mort_mat, 'year', 'timesteps')
mort_mat<- data.table(mort_mat)
mort_mat<- mort_mat[,timesteps:= timesteps* 365]
mort_mat <- mort_mat |>
  select(-timesteps)


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
intvn<- set_vaccine_coverage(intvn, 
                             change= TRUE,
                             terminal_year= 2050, 
                             rtss_target= 0.8, 
                             rtss_year= 2023)

baseline<- set_vaccine_coverage(baseline,
                                change= FALSE, 
                                terminal_year= 2050, 
                                rtss_target= 0.8, 
                                rtss_year= 2023)

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


# prep central burden estimate inputs ------------------------------------------

bl<- prep_inputs(baseline, 
                    mort_dt= mort, 
                    death_rate_matrix= mort_mat,
                    folder= paste0(malaria_dir, '/central_estimates/baseline/'))

int<- prep_inputs(intvn, 
                  mort_dt= mort, 
                  death_rate_matrix= mort_mat,
                  folder= paste0(malaria_dir, '/central_estimates/intervention/'))


# prep stochastic burden estimate inputs
int_stochastic<- lapply(int, prep_stochastic_inputs, draws= 10)
int_stochastic<- flatten(int_stochastic)

bl_stochastic<- lapply(bl, prep_stochastic_inputs, draws= 10)
bl_stochastic<- flatten(bl_stochastic)

# submit jobs to cluster  ------------------------------------------------------
message(paste0('submitting ', length(bl),  ' central burden jobs'))
message(paste0('submitting ', length(int),  ' central burden jobs'))


# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation")

# save a context (working environment for your code) ---------------------------
# additional script contains helper functions for larger scale model runs
ctx <- context::context_save('pkgs', 
                             packages = packages, 
                             package_sources = src,
                             sources = 'Q:/VIMC_malaria/run_malaria_model.R')

config <- didehpc::didehpc_config(use_rrq = TRUE)

# load context into queue
obj <- didehpc::queue_didehpc(ctx)
didehpc::web_login()

# run central burden baseline jobs ---------------------------------------------
central_fold<- paste0(malaria_dir, 
                      '/central_estimates/baseline/') # folder you would like to save outputs in
stochastic_fold<-  paste0(malaria_dir, 
                          '/stochastic_estimates/baseline/')

central_baseline_jobs_5k <- obj$lapply(bl_5k, 
                                    run_malaria_model, 
                                    folder= central_fold,
                                    stochastic_run= FALSE)

central_baseline_jobs_50k <- obj$lapply(bl_50k, 
                                       run_malaria_model, 
                                       folder= central_fold,
                                       stochastic_run= FALSE)
#stochastic_baseline_jobs<- obj$lapply(bl_stochastic, 
#run_malaria_model, 
#folder= stochastic_fold,
#stochastic_run= T)



# run central burden  intervention jobs ----------------------------------------
central_fold<- paste0(malaria_dir, 
                      '/central_estimates/intervention/') # folder you would like to save outputs in
stochastic_fold<-  paste0(malaria_dir, 
                          '/stochastic_estimates/intervention/')

central_intvn_jobs <- obj$lapply(int, 
                                 run_malaria_model, 
                                 folder= central_fold,
                                 stochastic_run= FALSE)

#stochastic_baseline_jobs<- obj$lapply(int_stochastic,  
#run_malaria_model, 
#folder= stochastic_fold)
test<- lapply(int,
              run_malaria_model, 
              folder= central_fold)

t <- obj$enqueue(1+5)
