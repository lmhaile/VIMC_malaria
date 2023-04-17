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
library(openxlsx)
library(wesanderson)
library(extrafont)

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/' #  directory where code is stored
malaria_dir<- 'F:/Lydia/VIMC'      #  project directory where files are stored
setwd('Q:/')

# load in inputs from VIMC  ----------------------------------------------------
pop<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_int_pop_both.csv'))
# le<- read.csv(paste0(malaria_dir, '/inputs/202212rfp-1_dds-202208_life_ex_both.csv'))
coverage<- read.csv(paste0(malaria_dir, '/inputs/coverage_202212rfp-1_malaria-mal4-default.csv'))
mort<- read.csv(paste0(malaria_dir, '/inputs/mortality.csv'))
# fert<- read.csv(paste0(code_dir, '/inputs/fertility.csv'))

# pull site data  --------------------------------------------------------------
mli<- foresite::MLI
site<- site::single_site(mli, 4)

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

# set vaccine coverage ------------------------------------------------------
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
                 death_rate_matrix= mort_mat)

int<- prep_inputs(intvn, 
                  mort_dt= mort, 
                  death_rate_matrix= mort_mat)

# prep stochastic burden estimate inputs
int_stochastic<- lapply(int, prep_stochastic_inputs, draws= 50)
bl_stochastic<- lapply(bl, prep_stochastic_inputs, draws= 50)


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
                             sources = 'Q:/VIMC/run_malaria_model.R')


# load context into queue
obj <- didehpc::queue_didehpc(ctx)

# run central burden baseline jobs ---------------------------------------------
central_fold<- paste0(malaria_dir, 
                      '/central_estimates/baseline/') # folder you would like to save outputs in
stochastic_fold<-  paste0(malaria_dir, 
                          '/stochastic_estimates/baseline/')

central_baseline_jobs <- obj$lapply(bl, 
                                    run_malaria_model, 
                                    folder= central_fold)

stochastic_baseline_jobs<- obj$lapply(bl_stochastic, 
                                      run_malaria_model, 
                                      folder= stochastic_fold)

# run central burden  intervention jobs ----------------------------------------
central_fold<- paste0(malaria_dir, 
                      '/VIMC/central_estimates/intervention/') # folder you would like to save outputs in
stochastic_fold<-  paste0(malaria_dir, 
                          '/VIMC/stochastic_estimates/intervention/')

central_baseline_jobs <- obj$lapply(int, 
                                    run_malaria_model, 
                                    folder= central_fold)

stochastic_baseline_jobs<- obj$lapply(int_stochastic, 
                                      run_malaria_model, 
                                      folder= stochastic_fold)


# load in files  ---------------------------------------------------------------
dir<- paste0(malaria_dir, 
             '/VIMC/central_estimates/baseline/') #directory where outputs are
files<- list.files(dir, full.names = T)
bl<- rbindlist(lapply(files, readRDS), fill= T)


dir<- paste0(malaria_dir, 
             '/VIMC/central_estimates/intervention/') #directory where outputs are
files<- list.files(dir, full.names = T)
intvn<- rbindlist(lapply(files, readRDS), fill= T)#

intvn<- intvn |>
  mutate(run = 'intervention')

bl<- bl |>
  mutate(run = 'baseline')

# aggregate model outputs  ----------------------------------------
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
write.csv(intvn, paste0(malaria_dir, '/VIMC/output/central_burden_estimates/central_burden_vaccine.csv'))
write.csv(bl, paste0(malaria_dir, '/VIMC/output/central_burden_estimates/central_burden_baseline.csv'))


# plot outputs over time  ------------------------------------------------------
intvn<- intvn |> mutate( scenario = 'intervention')
bl<- bl |> mutate(scenario = 'baseline')

output<- rbind(intvn, bl, fill= T)
#site_name<- unique(output$)
#font_import()
#loadfonts(device = 'win')

# clinical cases  --------------------------------------------------------------
ggplot(data= output, mapping = aes(x= year, y= cases, color= scenario, fill= scenario))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age) +
  labs(x= 'Time (in years)', y= 'Clinical cases', title= paste0('Clinical cases over time: ',site_name),
       color= 'Scenario', fill= 'Scenario') +
  theme_minimal()+
  theme(text= element_text(family= 'Calibri')) +
  scale_color_manual(values= wes_palette('Royal2', n= 2)) +
  scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
  
# deaths  ----------------------------------------------------------------------
ggplot(data= output, mapping = aes(x= year, y= deaths, color= scenario, fill= scenario))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age) +
  labs(x= 'Time (in years)', y= 'Deaths', title= 'Deaths over time',
       color= 'Scenario', fill= 'Scenario') +
  theme_minimal()+
  theme(text= element_text(family= 'Calibri')) +
  scale_color_manual(values= wes_palette('Royal2', n= 2)) +
  scale_fill_manual(values= wes_palette('Royal2', n= 2)) 


# DALYs ------------------------------------------------------------------------
ggplot(data= output, mapping = aes(x= year, y= dalys, color= scenario, fill= scenario))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age) +
  labs(x= 'Time (in years)', y= 'DALYs', title= 'DALYs over time',
       color= 'Scenario', fill= 'Scenario') +
  theme_minimal()+
  theme(text= element_text(family= 'Calibri')) +
  scale_color_manual(values= wes_palette('Royal2', n= 2)) +
  scale_fill_manual(values= wes_palette('Royal2', n= 2)) 


