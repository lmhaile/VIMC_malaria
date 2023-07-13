################################################################################
##  title   02_postprocess_outputs
##  author  Lydia Haile
##  purpose launches jobs to aggregate outputs, produce cases/deaths/dalys
################################################################################

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

# custom functions -------------------------------------------------------------
source("Q:/VIMC_malaria/functions/postprocessing_functions.R", echo=TRUE)

# directories ------------------------------------------------------------------

drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/'                   #  directory where code is stored
malaria_dir<- 'Q:/VIMC/central_estimates/'      #  project directory where files are stored
setwd('Q:/')


# load in files  ---------------------------------------------------------------
iso<- 'NGA'                                     # country you would like to process results for
tag<- 'run_through_2020'                      # description of the model run you are carrying out


# 
file_dir<- paste0(malaria_dir, 'raw_model_output/', iso, '/', tag, '/')
files<- list.files(file_dir, full.names= TRUE)[1]
input<- rbindlist(lapply(files, readRDS))

# dt <- get_rates(
#   dt,
#   time_divisor = 365,
#   baseline_t = 0,
#   age_divisor = 1,
#   scaler = 0.215,
#   treatment_scaler = 0.5,
#   baseline_treatment = 0
# )

#dt <- drop_burnin(input, burnin= 5*365)                         # drop the burn-in period from the data-set
dt <- time_transform(x= dt, time_divisor = 365, baseline_t = 0) # transform time into annual outputs
dt <- aggregate_outputs(dt, interval= 365)                      # aggregate model outputs  
dt <- calculate_deaths_ylls(dt)                                 # calculate deaths 
dt <- calculate_ylds_dalys(dt)                                  # calculate DALYs
dt <- reformat_vimc_outputs(dt)                                 # format outputs for submission


# check model output prevalence is similar to underlying prevalence for time period
# Calculate prevalence
dt<- copy(input)
dt$prevalence <- dt$n_detect_730_3649 / dt$n_730_3649


# Set the time
dt$t <- (dt$timestep / 365) + 2000

# Plot
plot(dt$prevalence ~ dt$t, t = "l", ylim = c(0, 0.8), xlab = "Year", ylab = "Prevalence")

# Add MAP prevalence
points(site_data$prevalence$year + 0.5, site_data$prevalence$pfpr, pch = 19, col = "darkred")

dev.off()

# save outputs
write.csv(intvn, paste0(malaria_dir, '/output/central_burden_estimates/central_burden_vaccine.csv'))
write.csv(bl, paste0(malaria_dir, '/output/central_burden_estimates/central_burden_baseline.csv'))
