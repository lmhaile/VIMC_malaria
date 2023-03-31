################################################################################
##  title   VIMC model test case
##  author  Lydia Haile
##  purpose carry out modelling done for VIMC RFP using malariasimulation package,
##          with one site to debug + troubleshoot
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

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/' #  directory where code is stored
malaria_dir<- 'F:/Lydia'      #  project directory where files are stored

# pull site data  --------------------------------------------------------------
mli<- foresite::MLI
site<- site::single_site(mli, 4)

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

# modify vaccine coverage -------------------------
terminal_year<- 2100 # year you would like to expand intervention coverage out to 

rtss_change<- T   # do you want to modify RTSS?
rtss_target<- 0.8 # target for RTSS
rtss_year<-  2023    # year for target

group_var <- names(site$sites)

# expand intervention years ----------------------------------------------------
intvn$interventions <- intvn$interventions |>
  expand_interventions(max_year = terminal_year,
                       group_var = group_var)

baseline$interventions <- baseline$interventions |>
  expand_interventions(max_year = terminal_year,
                       group_var = group_var)
# RTSS coverage ----------------------------------------------------------------
if (rtss_change== T){
  
  intvn$interventions <- intvn$interventions |>
    set_change_point(sites = intvn$sites, 
                     var = "rtss_cov", 
                     year = rtss_year, 
                     target = rtss_target)
}


# Linear scale up of coverage
intvn$interventions <- intvn$interventions |>
  linear_interpolate(vars = c("itn_use", "pmc_cov", "smc_cov", "rtss_cov"), 
                     group_var = group_var)

baseline$interventions <- baseline$interventions |>
  linear_interpolate(vars = c("itn_use", "pmc_cov", "smc_cov", "rtss_cov"), 
                     group_var = group_var)

intvn$interventions <- intvn$interventions |>
  fill_extrapolate(group_var = group_var)

baseline$interventions <- baseline$interventions |>
  fill_extrapolate(group_var = group_var)

baseline$interventions <- baseline$interventions |>
  add_future_net_dist(group_var = group_var)

intvn$interventions <- intvn$interventions |>
  add_future_net_dist(group_var = group_var)

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

# prep site data for model launch ----------------------------------------------
prep_inputs<- function(site_data){
  
  #' Prep inputs for batch launch
  #'
  #' @param site_data dataset with site files for country
  #' output: list with site name, urban/rural grouping, iso code, and parameters to pass into cluster
  
  
  # how many sites in this country?
  jobs<- nrow(site_data$sites)
  
  message(paste0('prepping ', jobs, ' jobs for model launch'))
  
  prep_site_data<- function(num){
    site<- site::single_site(site_file= site_data, index= num) 
    
    ## get site info
    site_name<- site$sites$name_1
    ur<- site$sites$urban_rural
    iso<- site$sites$iso3c
    message(paste0('prepping inputs for site ', site_name, ' ', ur))
    
    # pull parameters for this site
    params<- site::site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir= site$eir$eir[1],
      overrides = list(human_population= 10000) # what size population is appropriate?
    )
    
    inputs<- list('param_list'= params, 'site_name'= site_name, 'ur'= ur, 'iso'= iso)
    return(inputs)
  }
  output<- lapply(c(1:jobs), prep_site_data)
}


baseline<- prep_inputs(baseline)
intervention<- copy(intvn)


# submit jobs to cluster  ------------------------------------------------------
message(paste0('submitting ', length(output),  ' jobs'))

# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation")

# save a context (working environment for your code) ---------------------------
# additional script contains helper functions for larger scale model runs
ctx <- context::context_save('pkgs', 
                             packages = packages, 
                             package_sources = src,
                             sources = 'Q:/model_onboarding/run_malaria_model.R')


# load context into queue
obj <- didehpc::queue_didehpc(ctx)


# run baseline jobs
fold<- 'Q:/model_test_run/baseline/' # folder you would like to save outputs in
dir.create(fold)
grp1 <- obj$lapply(output, run_malaria_model, folder= fold)

# run intervention jobs
fold<- 'Q:/model_test_run/intervention/' # folder you would like to save outputs in
dir.create(fold)

grp2 <- obj$lapply(intervention, run_malaria_model, folder= fold)

