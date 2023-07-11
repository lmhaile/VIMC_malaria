# RFP comparison

rm(list= ls())

library(remotes)
library(didehpc)
library('drat')
library(readr)
library('tidyr')
library('data.table')
library(malariasimulation)
library(dplyr)
library(ggplot2)
library(ggforce)
library(readr)
source('Q:/VIMC_malaria/VIMC_functions.R')
source('Q:/VIMC_malaria/run_malaria_model.R')

dir<- 'Q:/VIMC_malaria_rfp/' #directory issue
save_dir<- 'Q:/VIMC_files/central_estimates/'

year <- 365
month <- 30

# Simulation length (50 years)
timesteps <- year * 50


# # run model on matched EIR ---------------------------------------------------
# 1.9 (PFPR 10) or 27.3 (PFPR 50)
p <- get_parameters(list(
  human_population = 500000,
  individual_mosquitoes = FALSE))


# Set initial level of transmission (change init EIR to match transmission to RFP runs)
# in RFP outputs EIR was closer to 3.12 over study period
p <- set_equilibrium(p, init_EIR = 27.3)

# Set up treatment
p <- set_drugs(p, list(AL_params))
p <- set_clinical_treatment(p, 1, 1, 0.5)


# Set up age-outputs
min_ages <- year * seq(0, 99, by= 1)
max_ages <- year * seq(1, 100, by= 1) -1

p$clinical_incidence_rendering_min_ages = min_ages
p$clinical_incidence_rendering_max_ages = max_ages
p$severe_incidence_rendering_min_ages = min_ages
p$severe_incidence_rendering_max_ages = max_ages

# Set age group rendering
p$age_group_rendering_min_ages = min_ages
p$age_group_rendering_max_ages = max_ages


# Set up treatment
p <- set_drugs(p, list(AL_params))
p <- set_clinical_treatment(p, 1, 1, 0.5)


# Set up vaccine for intvn_scenario
p_intvn <- set_pev_epi(
  p,
  profile = rtss_profile, # We will model implementation of the RTSS vaccine.
  timesteps = 35 * year, # Vaccination will begin at 20 years into the simulation.
  coverages = 0.8, # Vaccine coverage is 80%.
  min_wait = 0, # There is no minimum wait since the last vaccination.
  age = 6 * month, # Individuals will be vaccinated once they reach 6 months of age.
  booster_timestep = 18 * month, # The booster is administered 18 months following the third dose.
  booster_coverage = 0.8, # 80% of those vaccinated with the primary series will be boosted.
  booster_profile = list(rtss_booster_profile)
  )


# save input parameters to input folder
write_rds(p_intvn, paste0(save_dir, 'intervention/input_parameters_vaccine_RFP_PFPR_50.rds'))
write_rds(p, paste0(save_dir, 'baseline/input_parameters_baseline_RFP_PFPR_50.rds'))

# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation@dev")

# save a context (working environment for your code) ---------------------------
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
  template = "GeneralNodes",
  ## use for the wpia cluster
  parallel = FALSE
)

#load context into queue
obj <- didehpc::queue_didehpc(ctx, config)

#obj <- didehpc::queue_didehpc(ctx, provision = "later")
# obj$install_packages("mrc-ide/individual")
# obj$install_packages("mrc-ide/malariaEquilibrium")
# obj$install_packages("mrc-ide/malariasimulation@dev")


didehpc::web_login()


# run jobs ---------------------------------------------

# 10 % PFPR
bl_10 <-
  list(
    'parameter_filepath' = paste0(save_dir, 'baseline/input_parameters_baseline_RFP_PFPR_10.rds'),
    'output_folder' = 'Q:/VIMC_files/central_estimates/baseline/raw_model_output/',
    'tag' = 'PFPR_10'
  )

int_10 <-
  list(
    'parameter_filepath' =  paste0(save_dir, 'intervention/input_parameters_vaccine_RFP_PFPR_10.rds'),
    'output_folder' = 'Q:/VIMC_files/central_estimates/intervention/raw_model_output/',
    'tag' = 'PFPR_10'
  )
# 50 % PFPR
bl_50 <-
  list(
    'parameter_filepath' = paste0(save_dir, 'baseline/input_parameters_baseline_RFP_PFPR_50.rds'),
    'output_folder' = 'Q:/VIMC_files/central_estimates/baseline/raw_model_output/',
    'tag'= 'PFPR_50'
  )

int_50 <-
  list(
    'parameter_filepath' =  paste0(save_dir, 'intervention/input_parameters_vaccine_RFP_PFPR_50.rds'),
    'output_folder' = 'Q:/VIMC_files/central_estimates/intervention/raw_model_output/',
    'tag' = 'PFPR_50'
  )


filled<- list(bl_10, int_10, bl_50, int_50)

test <- obj$lapply(filled, run_malaria_model_rfp)

# test locally  ---------------------------------------------------------------
run_malaria_model_rfp(input_filepath= paste0(save_dir, 'baseline/input_parameters_baseline_RFP.rds'),
                      folder= 'Q:/VIMC_files/central_estimates/baseline/raw_model_output/')

run_malaria_model_rfp(input_filepath= paste0(save_dir, 'intervention/input_parameters_vaccine_RFP.rds'),
                      folder= 'Q:/VIMC_files/central_estimates/intervention/raw_model_output/')


write.xlsx(model, 'Q:/VIMC_files/central_estimates/intervention/raw_model_output/rfp_output.xlsx')


