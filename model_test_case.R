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

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/' #  directory where code is stored
malaria_dir<- 'F:/Lydia'      #  project directory where files are stored
setwd('Q:/')

# load in inputs from VIMC  ----------------------------------------------------
pop<- read.csv(paste0(code_dir, 'inputs/202212rfp-1_dds-202208_int_pop_both.csv'))
le<- read.csv(paste0(code_dir, 'inputs/202212rfp-1_dds-202208_life_ex_both.csv'))
coverage<- read.csv(paste0(code_dir, 'inputs/coverage_202212rfp-1_malaria-mal4-default.csv'))
mort<- read.csv(paste0(code_dir, 'inputs/mortality.csv'))
fert<- read.csv(paste0(code_dir, 'inputs/fertility.csv'))

# pull site data  --------------------------------------------------------------
mli<- foresite::MLI
site<- site::single_site(mli, 4)

# format mortality data
mort<- mort[, c('age_to', 'year', 'value')]
mort<- data.table(mort)
mort[age_to== 0, age_to:= 1]
mort[, age_to:= age_to * 365]

#bind on a baseline year
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

# modify vaccine coverage -------------------------
terminal_year<- 2050 # year you would like to expand intervention coverage out to 
# we would like to run model out to 2100, but do not have population projections for 2050-2100 yet

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
prep_inputs<- function(site_data, mort_dt, death_rate_matrix){
  
  #' Prep inputs for batch launch for GAVI runs
  #'
  #' @param site_data dataset with site files for country
  #' @param mort_dt   
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
      min_ages = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95) *365,
      eir= site$eir$eir[1],
      overrides = list(human_population= 10000) # what size population is appropriate?
    )
    
    # # set custom demography based on mortality inputs ------------------------
    # params<- set_demography(
    #   params,
    #   agegroups= unique(mort_dt$age_to),
    #   timesteps = unique(mort_dt$year)*365,
    #   deathrates = death_rate_matrix
    #   
    # )
    # 
    inputs<- list('param_list'= params, 'site_name'= site_name, 'ur'= ur, 'iso'= iso)
    return(inputs)
  }
  output<- lapply(c(1:jobs), prep_site_data)
}


bl<- prep_inputs(baseline, mort_dt= mort, death_rate_matrix= mort_mat)

int<- prep_inputs(intvn, mort_dt= mort, death_rate_matrix= mort_mat)


# submit jobs to cluster  ------------------------------------------------------
message(paste0('submitting ', length(bl),  ' jobs'))
message(paste0('submitting ', length(int),  ' jobs'))

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

# run baseline jobs
fold<- paste0(malaria_dir, '/VIMC/baseline/') # folder you would like to save outputs in
dir.create(fold)
grp1 <- obj$lapply(bl, run_malaria_model, folder= fold)

# run intervention jobs
fold<- paste0(malaria_dir, '/VIMC/intervention/') # folder you would like to save outputs in
dir.create(fold)
grp2 <- obj$lapply(int, run_malaria_model, folder= fold)

lapply(bl, run_malaria_model, folder= fold)


# load in files  ------------------------------------------
dir<- paste0(malaria_dir, '/VIMC/baseline/') #directory where outputs are
files<- list.files(dir, full.names = T)
bl<- rbindlist(lapply(files, readRDS), fill= T)


dir<- paste0(malaria_dir, '/VIMC/intervention/') #directory where outputs are
files<- list.files(dir, full.names = T)
intvn<- rbindlist(lapply(files, readRDS), fill= T)#

intvn<- intvn |>
  mutate(run = 'intervention')

bl<- bl |>
  mutate(run = 'baseline')

# reformat and aggregate model outputs  ----------------------------------------

aggregate_outputs<- function(dt, interval){
  
  #' aggregate incident cases based on a pre-determined time interval (expressed in days). 
  #' sum to the country level.
  #' 
  #' @param dt       raw model output from malariasimulation package
  #' @param interval time period you would like to calculate incidence over, expressed in days.
  #' 
  #' output: data table with summed cases and rates over specified time interval.
  
  # reformat case outputs to long
  # need clinical cases, severe cases, population, number treated, and number detected
  
  message(paste0('aggregating outputs by time interval: ', interval, ' days'))
  dt <- dt |> 
    select(timestep, iso, run,
           contains("n_inc_clin"), contains("n_inc_sev"), contains("n_age"), contains('n_treated'), contains('n_detect')) |>  
    pivot_longer(c(contains("n_inc_clin"), contains("n_inc_sev"), contains("n_age")),
                 names_to = "age", 
                 values_to = "value") |>
    mutate(type = ifelse(grepl('n_inc_clinical', age), 'clinical', 'severe')) |>
    mutate(type = ifelse(grepl('n_age', age), 'population', type))|>
    #mutate(type = ifelse(grepl('n_detect', age), 'detected', type))|>
    mutate(age = gsub('n_inc_clinical_', '', age),
           age = gsub('n_inc_severe_', '', age),
           age = gsub('n_age_', '', age)) |>
    spread(key= type, value= value) |>
    separate(age, into = c("age_days_start", "age_days_end"), sep = "_") 
  
  
  # calculate incidence based on some time interval
  dt <- dt |>
    mutate(time = as.integer(timestep/ interval))
  
  
  # sum cases based on this interval, also by country
  dt<- dt |> 
    group_by(age_days_start, time, iso) |> 
    mutate(clinical = sum(clinical),
           severe = sum(severe),
           n_treated = sum(n_treated),
           #n_detect= sum(detected),
           population = round(mean(population))) |>
    select(-timestep, -n_detect_365_36499, -n_detect_730_3649) |>
    distinct()
  
  
  
  # calculate rates based on this interval
  dt<- dt |> 
    mutate(clin_rate = clinical/ population,
           severe_rate = severe/ population)
  #prevalence= n_detect/ population)
  
  return(dt)
  message('completed aggregation')
  
}

intvn<-aggregate_outputs(intvn, interval= 30)
bl<-aggregate_outputs(bl, interval= 30)


# calculate deaths -------------------------------------------------------------
# transform age into years

calculate_deaths_ylls<- function(dt, cfr= 0.215, treatment_scaler= 0.5, lifespan= 63){
  
  #' Calculate deaths + years of life lost (YLLs) per GTS method.
  #' Where severe cases= 0, deaths= 0. Additionally remove a proportion of the cases that have received treatment.
  #' 
  #' @param dt               malariasimulation model outputs with columns 'severe' for severe incidence and 'n_treated' for number of individuals treated
  #' @param cfr              per GTS method, deaths are calculated using a case fatality ratio (CFR) value that is applied to severe incidence.
  #'                         see World Malaria Report and Wilson et al. for more information.
  #' @param treatment_scaler we remove a proportion of cases that have received treatment, assuming that 50% of treated cases remit and
  #'                         are no longer susceptible to mortality.
  #' @param lifespan         expected lifespan used to calculate YLLs. 
  #'                         YLLs are calculated by multiplying deaths by the number of years an individual was expected to live past their year of death.
  #'                         when comparing YLLs across different locations, it is recommended to use the same lifespan across YLL calculations.
  #'                         Typically, you should use the highest observed life expectancy in the region/ location you are studying.
  #' Output: data table with columns titled 'deaths' and 'yll'
  
  
  message('calculating deaths and YLLs')
  
  # transform age into years
  dt<- dt|> 
    mutate(age_days_start= as.numeric(age_days_start),
           age_days_end= as.numeric(age_days_end))
  
  dt<- dt |>
    mutate(age_years_start= as.integer(age_days_start/ 365),
           age_years_end= as.integer(age_days_end/ 365))
  dt<- dt |>
    mutate(deaths =  cfr * severe - (n_treated * treatment_scaler))
  
  
  dt<- data.table(dt)
  
  dt[deaths< 0 , deaths:= 0]
  dt[, deaths:= as.integer(deaths)]
  
  # calculate ylls
  
  dt <- dt |>
    mutate(yll= deaths * (lifespan - (age_years_end - age_years_start)/2))  
  
  message('completed calculation of deaths and YLLs')
  
  return(dt)
}


intvn<- calculate_deaths_ylls(intvn)
bl<- calculate_deaths_ylls(bl)

calculate_ylds_dalys<- function(dt, 
                                mild_dw= 0.006, 
                                moderate_dw= 0.051, 
                                severe_dw= 0.133,
                                clin_episode_length= 0.01375,
                                severe_episode_length= 0.04795){
  
  #' Calculate Years Lived with Disability (YLDs) and Disability-Adjusted Life-Years
  #' based on disability weights from the Global Burden of Disease study.
  #' 
  #' Disability weights sourced here:
  #' https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fghdx.healthdata.org%2Fsites%2Fdefault%2Ffiles%2Frecord-attached-files%2FIHME_GBD_2017_DISABILITY_WEIGHTS_Y2018M11D08.XLSX&wdOrigin=BROWSELINK
  #' 
  #' Keep in mind this is an approximation of YLD estimation from the GBD study; disability due to comorbid conditions
  #' such as motor impairment and anemia are excluded.
  #' 
  #' For now, we assume that cases under 5 are moderate, due to the higher severity of malaria at younger ages.
  #' This assumption may be revisited in the future, 
  #' @param dt                    input dataset with columns 'clinical' for clinical incidence,
  #'                              'severe' for severe incidence, 'deaths', 'ylls', 'age_years_start', and 'age_years_end'
  #' @param mild_dw               disability weight for mild malaria
  #' @param moderate_dw           disability weight for moderate malaria
  #' @param severe_dw             disability weight for severe malaria
  #' @param clin_episode_length   length of an episode of clinical malaria
  #' @param severe_episode_length length of an episode of severe malaria 
  #' 
  #' Output: data table with columns 'ylds' and 'dalys'
  
  require(data.table)
  
  message('calculating YLDs and DALYs')
  
  # calculate YLDs  ---
  dt[ age_years_end < 5, 
      yld:= severe * severe_dw * severe_episode_length + 
        clinical * moderate_dw * clin_episode_length]
  
  dt[ age_years_end >= 5, 
      yld:= severe * severe_dw * severe_episode_length + 
        clinical * mild_dw * clin_episode_length]
  
  # calculate DALYs ---
  dt<- dt |>
    mutate(daly= yll+ yld)
  
  message('completed calculation of YLDs and DALYs')
  
  return(dt)
}


intvn<- calculate_ylds_dalys(intvn)
bl<- calculate_ylds_dalys(bl)


calculate_cases_deaths_averted<- function(intvn_dt, bl_dt){
  #' Calculates cases, deaths, and dalys averted between intervention and baseline scenario.
  #' @param intvn_dt               model outputs for intervention scenario
  #' @param baseline               model outputs for baseline scenario
  #' 

  intvn_dt<- intvn_dt |>
    select(time, iso, run, age_years_start, age_years_end,
           clinical, severe, deaths, daly) |>
    rename(clinical_intvn= clinical,
           severe_intvn= severe,
           deaths_intvn= deaths,
           daly_intvn= daly)
  
  bl_dt<- bl_dt |>
    select(time, iso, run, age_years_start, age_years_end,
           clinical, severe, deaths, daly) |>
    rename(clinical_bl= clinical,
           severe_bl= severe,
           deaths_bl= deaths,
           daly_bl= daly)
  
  output<- left_join(bl_dt, intvn_dt, 
                     by= c('time', 'iso', 'age_years_start', 'age_years_end'))
  
  # calculate averted cases/deaths/dalys
  
  output<- output|>
    mutate(cases_averted = clinical_bl - clinical_intvn,
           severe_cases_averted = severe_bl - severe_intvn,
           deaths_averted = deaths_bl - deaths_intvn,
           dalys_averted = daly_bl - daly_intvn)
  
  message('calculated cases deaths dalys averted')
  
  return(output)
  
}
# summarize
# dt |> 
#  group_by(age_years_start) |> 
#  summarise(yll= sum(yll),
#            yld= sum(yld),
#            daly= sum(daly))
#

# plot outputs over time  ------------------------------------------------------
output<- rbind(intvn, bl, fill= T)

# clinical cases
ggplot(data= output, mapping = aes(x= time, y= clinical, color= run, fill= run))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age_years_start) +
  labs(x= 'Time (in months)', y= 'Clinical cases (by months)', main= 'Clinical cases over time') +
  theme_minimal()

# severe cases
ggplot(data= output, mapping = aes(x= time, y= severe, color= run, fill= run))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age_years_start) +
  labs(x= 'Time (in months)', y= 'Severe cases (by months)', main= 'Severe cases over time') +
  theme_minimal()

# deaths
ggplot(data= output, mapping = aes(x= time, y= deaths, color= run, fill= run))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age_years_start) +
  labs(x= 'Time (in months)', y= 'Deaths (by months)', main= 'Deaths over time') +
  theme_minimal()

# dalys
ggplot(data= output, mapping = aes(x= time, y= daly, color= run, fill= run))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age_years_start) +
  labs(x= 'Time (in months)', y= 'DALYs (by months)', main= 'DALYs over time') +
  theme_minimal()

