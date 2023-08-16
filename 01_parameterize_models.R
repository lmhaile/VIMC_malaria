################################################################################
##  title   01_model_launch
##  author  Lydia Haile
##  purpose parameterize models and launch on HPC cluster
################################################################################

rm(list= ls())
setwd('C:/Users/lhaile/Documents/VIMC_malaria')

# packages  --------------------------------------------------------------------
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
install.packages("foresite", type = "source") # v0.1.0
# custom functions to source  --------------------------------------------------
source("functions/parameterizing_functions.R")
source("functions/modelling_functions.R")

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/'      #  directory where code is stored
malaria_dir<- 'Q:/VIMC/central_estimates/'      #  project directory where files are stored
input_dir<- 'Q:/VIMC'
setwd('Q:/')

# load in inputs from VIMC  ----------------------------------------------------
coverage<- read.csv(paste0(input_dir, '/inputs/coverage_202212rfp-1_malaria-mal4-default.csv'))
mort<- read.csv(paste0(input_dir, '/inputs/mortality.csv'))
countries <- read.csv("M:/Lydia/malaria-sites/raw_data/wmr/wmr_countries.csv")
iso<- countries$iso3c

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

# plot initial intervention coverage  -----------------------------------------------
site_data<- foresite::NGA

year<- 365
iso<- 'NGA'
descrip<- 'fixed_demography_test_run'
cluster<- TRUE #launch on cluster? 

# for large scale runs of full countries  --------------------------------------

# launch parameterization on cluster
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation', 'foresite')
src <- conan::conan_sources("github::mrc-ide/malariasimulation", repos = "https://mrc-ide.github.io/drat/")

ctx <- context::context_save(
  'pkgs',
  packages = packages,
  package_sources = src,
  sources = 'functions/parameterizing_functions.R'
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

jobs<- obj$lapply(filepaths, run_malaria_model)

  prep_country(iso,
               scenario= 'intervention', 
               vimc_mortality= TRUE,
               mortality_input = mort,
               min_ages = c(seq(0, 14, by= 1), seq(15, 80, by= 15)) * year,
               max_ages = c(seq(1, 15, by= 1), seq(30, 95, by= 15)) * year -1,
               population = 50000,
               description= descrip)




# plot demography


# for single site test runs ----------------------------------------------------
site_data<- foresite::NGA
site_data<- site::single_site(site_data, 2)

prep_site(site_data,
          scenario= 'intervention', 
          vimc_mortality= TRUE,
          mortality_input = mort,
          min_ages = c(seq(0, 14, by= 1), seq(15, 80, by= 15)) * year,
          max_ages = c(seq(1, 15, by= 1), seq(30, 95, by= 15)) * year -1,
          population = 10000,
          description= descrip)


demography_diagnostic(site_data)
demography_diagnostic(site)
# submit jobs to cluster -------------------------------------------------------
filepaths<- list.files(paste0(malaria_dir, 'input_parameters/', iso, '/', descrip),
                       full.names= T)

filepaths<- filepaths[1:10]
if (cluster== TRUE){
  packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
  src <- conan::conan_sources("github::mrc-ide/malariasimulation")
  
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
  
  jobs<- obj$lapply(filepaths, run_malaria_model)
  
}

# or launch a model locally  ---------------------------------------------------
filepath<- filepaths[[1]]
lapply(filepath, run_malaria_model)





# diagnostics ------------------------------------------------------------------
plot_interventions_combined(
  interventions = site_data$interventions,
  population = site_data$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment","SMC", "PMC")
)



