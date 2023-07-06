# RFP comparison

rm(list= ls())

#remotes::install_github("mrc-ide/malariasimulation@dev")
# remotes::install_github("mrc-ide/postie@ages")

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

human_population_val<- 5000

# read in input parameters for RFP
mwi_inputs<- readRDS(paste0(dir, "gavi_sites_inputs_draw0_mwi.rds"))

# Simulation length (50 years)
timesteps <- year * 50


# # first match PfPR to EIR for model launch  --------------------------------------------------------------
# # Set defualt parameters
# p <- get_parameters(list(
#   human_population = human_population_val,
#   prevalence_rendering_min_ages = 2 * year,
#   prevalence_rendering_max_ages = 10 * year,
#   individual_mosquitoes = FALSE))
# 
# 
# # Set up treatment
# p <- set_drugs(p, list(AL_params))
# p <- set_clinical_treatment(p, 1, 1, 0.5)
# 
# 
# # Set up age-outputs
# min_ages <- c(seq(0, 80, by = 20)) * year
# max_ages <- c(seq(20, 100, by = 20)) * year -1
# 
# p$clinical_incidence_rendering_min_ages = min_ages
# p$clinical_incidence_rendering_max_ages = max_ages
# p$severe_incidence_rendering_min_ages = min_ages
# p$severe_incidence_rendering_max_ages = max_ages
# 
# # Set age group rendering
# p$age_group_rendering_min_ages = min_ages
# p$age_group_rendering_max_ages = max_ages
# 
# 
# # Set up vaccine for intvn_scenario
# p_intvn <- set_pev_epi(
#   p,
#   profile = rtss_profile, # We will model implementation of the RTSS vaccine.
#   timesteps = 35 * year, # Vaccination will begin at 20 years into the simulation.
#   coverages = 0.8, # Vaccine coverage is 80%.
#   min_wait = 0, # There is no minimum wait since the last vaccination.
#   age = 6 * month, # Individuals will be vaccinated once they reach 6 months of age.
#   booster_timestep = 18 * month, # The booster is administered 12 months following the third dose.
#   booster_coverage = 0.8, # 95% of those vaccinated with the primary series will be boosted.
#   booster_profile = list(rtss_booster_profile))
# 
# 
# # Use the set_species() function to specify the mosquito population (species and
# # relative abundances)
# simparams <- set_species(parameters = p_intvn,
#                          species = list(arab_params, fun_params, gamb_params),
#                          proportions = c(0.25, 0.25, 0.5))
# 
# # Use the set_drugs() function to append the in-built parameters for the
# # drug Artemether Lumefantrine (AL)
# simparams <- set_drugs(simparams, list(AL_params))
# 
# # Use the set_clinical_treatment() function to parameterise human
# # population treatment with AL in the first timestep
# simparams <- set_clinical_treatment(parameters = simparams,
#                                     drug = 1,
#                                     timesteps = c(1),
#                                     coverages = c(0.45))
# 
# # Specify the time frame over which to simulate and the human population size:
# year <- 365
# human_population <- 5000
# 
# # Establish a vector of initial EIR values to simulate over and generate matching
# # PfPR2-10 values:
# init_EIR <- c(0.01, 0.1, 1, 5, 10, 25, 50, 60, 65, 75)
# 
# # For each initial EIR, calculate equilibrial parameter set and run the simulation
# malSim_outs <- lapply(
#   init_EIR,
#   function(init) {
#     p_i <- set_equilibrium(simparams, init)
#     run_simulation(5 * year, p_i)
#   }
# )
# 
# # Convert the default EIR output (per vector species, per timestep, across
# # the entire human population) to a cross-vector species average EIR per
# # person per year across the final year of the simulation:
# malSim_EIR <- lapply(
#   malSim_outs,
#   function(output) {
#     mean(
#       rowSums(
#         output[
#           output$timestep %in% seq(4 * 365, 5 * 365),
#           grepl('EIR_', names(output))
#         ] / human_population * year
#       )
#     )
#   }
# )
# # Calculate the average PfPR2-10 value across the final year for each initial
# # EIR value
# malSim_prev <- lapply(malSim_outs,
#                       function(output) {
#                         mean(output[output$timestep %in% seq(4 * 365, 5 * 365),
#                                     'n_detect_730_3650'] / output[output$timestep %in% seq(4 * 365, 5 * 365),
#                                                                   'n_730_3650'])
#                       })
# 
# # Create dataframe of initial EIR, output EIR, and prev 2-10 results
# malSim_P2E <-
#   cbind.data.frame(init_EIR,
#                    EIR = unlist(malSim_EIR),
#                    prev = unlist(malSim_prev))
# 
# 
# # Fit a line of best fit through malariasimulation initial EIR and prevalence
# malSim_fit <- predict(gam(prev ~ s(init_EIR), data = malSim_P2E),
#                       newdata = data.frame(init_EIR = c(0, seq(0.1, 75, 0.1))),
#                       type = "response")
# 
# # Append vector of initial EIR values to model fit to establish dataframe
# # for plotting
# malSim_fit <-
#   cbind(malSim_fit, data.frame(init_EIR = c(0 , seq(0.1, 75, 0.1))))
# 
# 
# # Define a colour palette for plotting:
# plot_cols <-
#   c("#E69F00",
#              "#56B4E9",
#              "#009E73",
#              "#CC79A7",
#              "#F0E442",
#              "#0072B2",
#              "#D55E00")
# 
# # Establish a plotting window:
# plot(x = 1, type = "n",
#      frame = F,
#      xlab = "Initial EIR", ylab = expression(paste(italic(Pf),"PR"[2-10])),
#      xlim = c(0,80), ylim = c(0, 1),
#      xaxs = "i", yaxs = "i")
# 
# # Overlay the initial EIR and corresponding PfPR2-10 points from malariasimulation
# points(x = malSim_P2E$init_EIR,
#        y = malSim_P2E$prev,
#        pch = 19,
#        col = 'black')
# 
# # Overlay the malariasimulation line of best fit:
# lines(x = malSim_fit$init_EIR,
#       y = malSim_fit$malSim_fit,
#       col = plot_cols[1],
#       lwd = 2,
#       type = "l",
#       lty = 1)
# 
# 
# # Store some pre-intervention baseline PfPR2-10 values
# PfPRs_to_match <- c(.10, .17, .35, .45, .50)
# 
# # Create a function to these baseline PfPR2-10 values to EIR values using
# # our model fit:
# match_EIR_to_PfPR <- function(x){
# 
#   m <- which.min(abs(malSim_fit$malSim_fit-x))
#   malSim_fit[m,2]
# 
# }
# 
# # Use the function to extract the EIR values for the specified
# # PfPR2-10 values:
# matched_EIRs <- unlist(lapply(PfPRs_to_match, match_EIR_to_PfPR))
# 
# # # # Create a dataframe of matched PfPR2-10 and EIR values:
# cbind.data.frame(PfPR = PfPRs_to_match, Matched_EIR = matched_EIRs)
# # # #  PfPR Matched_EIR
# # # #1 0.10         2.0
# # # #2 0.17         4.0
# # # #3 0.35        10.9
# # # #4 0.45        18.2
# # 
# 
# 
# # PfPR Matched_EIR
# # 1 0.10         1.9
# # 2 0.17         4.3
# # 3 0.35        12.2
# # 4 0.45        20.2
# # 5 0.50        27.3

# # using malariaequilibrium approach: EIR 1.6
# # run model on matched EIR ---------------------------------------------------
# 1.9 (PFPR 10) or 27.3 (PFPR 50)
p <- get_parameters(list(
  human_population = 50000,
  individual_mosquitoes = FALSE))


# Set initial level of transmission (change init EIR to match transmission to RFP runs)
# in RFP outputs EIR was closer to 3.12 over study period
p <- set_equilibrium(p, init_EIR = 1.9)

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
write_rds(p_intvn, paste0(save_dir, 'intervention/input_parameters_vaccine_RFP_PFPR_10.rds'))
write_rds(p, paste0(save_dir, 'baseline/input_parameters_baseline_RFP_PFPR_10.rds'))

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
filled<- list(bl_10, int_10)

test <- obj$lapply(filled, run_malaria_model_rfp)

# test locally  ---------------------------------------------------------------
run_malaria_model_rfp(input_filepath= paste0(save_dir, 'baseline/input_parameters_baseline_RFP.rds'),
                      folder= 'Q:/VIMC_files/central_estimates/baseline/raw_model_output/')

run_malaria_model_rfp(input_filepath= paste0(save_dir, 'intervention/input_parameters_vaccine_RFP.rds'),
                      folder= 'Q:/VIMC_files/central_estimates/intervention/raw_model_output/')


write.xlsx(model, 'Q:/VIMC_files/central_estimates/intervention/raw_model_output/rfp_output.xlsx')


