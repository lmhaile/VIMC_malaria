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
source("Q:/VIMC_malaria/VIMC_functions.R", echo=TRUE)

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/' #  directory where code is stored
malaria_dir<- 'Q:/VIMC_files'      #  project directory where files are stored
setwd('Q:/')


# load in files  ---------------------------------------------------------------
dir<- paste0(malaria_dir, 
             '/central_estimates/baseline/') #directory where outputs are
files<- list.files(dir, full.names = T)[1]
bl<- rbindlist(lapply(files, readRDS), fill= T)

dir<- paste0(malaria_dir, 
             '/central_estimates/intervention/') #directory where outputs are
files<- list.files(dir, full.names = T)
intvn<- rbindlist(lapply(files, readRDS), fill= T)

intvn<- intvn |>
  mutate(run = 'intervention')

bl<- bl |>
  mutate(run = 'baseline')

# drop the burn-in period from the data-set-- start with five years for now
intvn<- drop_burnin(intvn, burnin= 5*365)
bl<- drop_burnin(bl, burnin= 5*365)


# transform time into annual outputs
intvn<- time_transform(x= intvn, time_divisor = 365, baseline_t = 0)
bl<- time_transform(x= bl, time_divisor = 365, baseline_t = 0)


# aggregate model outputs  -----------------------------------------------------
# must be done at the site level first

sites<- unique(list(intvn$site_name))

test<- split(intvn, f = intvn$site_name)

testing<- lapply(test, aggregate_outputs, interval= 365, sum_to_country= F)


intvn<-aggregate_outputs(intvn, interval= 365)
bl<-aggregate_outputs(bl, interval= 365)


# calculate deaths -------------------------------------------------------------
intvn<- calculate_deaths_ylls(intvn)
bl<- calculate_deaths_ylls(bl)


# calculate DALYs --------------------------------------------------------------
intvn<- calculate_ylds_dalys(intvn)
bl<- calculate_ylds_dalys(bl)


# format outputs  --------------------------------------------------------------
intvn<- reformat_vimc_outputs(intvn)
bl<- reformat_vimc_outputs(bl)


# save output file to submission folder
write.csv(intvn, paste0(malaria_dir, '/output/central_burden_estimates/central_burden_vaccine.csv'))
write.csv(bl, paste0(malaria_dir, '/output/central_burden_estimates/central_burden_baseline.csv'))
