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
library(ggplot2)
library(ggthemes)
library(ggforce)
#install.packages('Q:/postie_0.1.2.tar.gz')
library(postie)


# custom functions -------------------------------------------------------------
#source("Q:/VIMC_malaria/functions/postprocessing_functions.R", echo=TRUE)

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/'                   #  directory where code is stored
malaria_dir<- 'Q:/VIMC/central_estimates/'      #  project directory where files are stored
setwd('Q:/')


# load in files  ---------------------------------------------------------------
iso<- 'NGA'
descrip<- 'no_changes_in_settings'

# make directories to save outputs to
raw_dir<- paste0(malaria_dir, 'raw_model_output/',  iso, '/', descrip, '/')
output_dir<- paste0(malaria_dir, 'processed_output/',  iso, '/', descrip, '/')

if(dir.exists(output_dir)== FALSE){
  dir.create(output_dir, recursive = T)
}

# directory where output files are  --------------------------------------------
files<- list.files(raw_dir, full.names= TRUE)

# generate cases, deaths, and dalys at the site level
output <- rbindlist(lapply(files, generate_vimc_output))                              
#write_rds(output, file= paste0(output_dir, 'site_level_output.RDS'))

# sum cases, deaths, and dalys to the country level
# aggregate<- aggregate_output(output)
# write_rds(aggregate, file= paste0(output_dir, 'country_level_output.RDS'))

# format for final submission
# final_output<- format_for_submission(aggregate)
# writeRDS(aggregate, file= paste0(output_dir, 'final_output.RDS'))

# diagnostics for a single site ------------------------------------------------
#plot_model_against_prevalence(input)

population_diagnostic(output)
incident_cases_diagnostic(output)
incidence_rate_diagnostic(output)
mortality_diagnostic(output)
mortality_rate_diagnostic(output)
daly_diagnostic(output)


# diagnostics for a single site  --------------------------------------------
population_diagnostic<- function(dt){
  
  p<-   ggplot(data= dt, mapping = aes(x= year, y= cohort_size, color= scenario, fill= scenario))+
    geom_point(alpha= 0.5)  +
    facet_wrap_paginate(~age) +
    labs(x= 'Time (in years)', y= 'Population', title= paste0('Population over time: site ', unique(dt$site_name)),
         color= 'Scenario', fill= 'Scenario') +
    theme_minimal()+
    theme(text= element_text(family= 'Arial Narrow')) +
    scale_color_manual(values= wes_palette('Royal2', n= 2)) +
    scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
  
  return(p)
  
}

incident_cases_diagnostic<- function(dt){
  
  p<-   ggplot(data= dt, mapping = aes(x= year, y= cases, color= scenario, fill= scenario))+
    geom_point(alpha= 0.5)  +
    facet_wrap_paginate(~age, scales= 'free') +
    labs(x= 'Time (in years)', y= 'Clinical cases', title= paste0('Incident clinical cases over time: site ', unique(dt$site_name)),
         color= 'Scenario', fill= 'Scenario') +
    theme_minimal()+
    theme(text= element_text(family= 'Arial Narrow')) +
    scale_color_manual(values= wes_palette('Royal2', n= 2)) +
    scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
  
  return(p)
  
}

incidence_rate_diagnostic<- function(dt){
  
  p<-  ggplot(data= dt, mapping = aes(x= year, y= clinical, color= scenario, fill= scenario))+
    geom_point(alpha= 0.5)  +
    facet_wrap_paginate(~age, 
                        scales = 'free') +
    labs(x= 'Time (in years)', y= 'Incidence rate', title= paste0('Incidence rate over time: ', unique(dt$site_name)),
         color= 'Scenario', fill= 'Scenario') +
    theme_minimal()+
    theme(text= element_text(family= 'Arial')) +
    scale_color_manual(values= wes_palette('Royal2', n= 2)) +
    scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
  
  return(p)
}



mortality_diagnostic<- function(dt){
  
  p<- ggplot(data= dt, mapping = aes(x= year, y= deaths, color= scenario, fill= scenario))+
    geom_point(alpha= 0.5)  +
    facet_wrap_paginate(~age, 
                        scales = 'free') +
    labs(x= 'Time (in years)', y= 'Deaths', title= paste0('Deaths over time: ', unique(output$site_name)),
         color= 'Scenario', fill= 'Scenario') +
    theme_minimal()+
    theme(text= element_text(family= 'Arial')) +
    scale_color_manual(values= wes_palette('Royal2', n= 2)) +
    scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
  
  return(p)
  
}
mortality_rate_diagnostic<- function(dt){
  
  p<- ggplot(data= dt, mapping = aes(x= year, y= mortality, color= scenario, fill= scenario))+
    geom_point(alpha= 0.5)  +
    facet_wrap_paginate(~age, 
                        scales = 'free') +
    labs(x= 'Time (in years)', y= 'Mortality rate', title= paste0('Mortality rate over time: ', unique(output$site_name)),
         color= 'Scenario', fill= 'Scenario') +
    theme_minimal()+
    theme(text= element_text(family= 'Arial')) +
    scale_color_manual(values= wes_palette('Royal2', n= 2)) +
    scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
  
  return(p)
  
}


daly_diagnostic<- function(dt){
  
  p<- ggplot(data= dt, mapping = aes(x= year, y= dalys, color= scenario, fill= scenario))+
    geom_point(alpha= 0.5)  +
    facet_wrap(~age, 
               scales = 'free') +
    labs(x= 'Time (in years)', y= 'DALYs', title= paste0('DALYs over time: ', unique(dt$site_name)),
         color= 'Scenario', fill= 'Scenario') +
    theme_minimal()+
    theme(text= element_text(family= 'Arial')) +
    scale_color_manual(values= wes_palette('Royal2', n= 2)) +
    scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
  
  return(p)
}



cases_averted_diagnostic<- function(dt){
  
  
  
}

# save outputs
#write.csv(intvn, paste0(malaria_dir, '/output/central_burden_estimates/central_burden_vaccine.csv'))
#write.csv(bl, paste0(malaria_dir, '/output/central_burden_estimates/central_burden_baseline.csv'))

