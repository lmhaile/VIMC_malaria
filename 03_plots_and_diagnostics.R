################################################################################
##  title   03_plots_and_diagnostics
##  author  Lydia Haile
##  purpose plot outputs
################################################################################

rm(list= ls())

# packages  --------------------------------------------------------------------
library(furrr)
library(data.table)
library(drat)
library(foresite)
library(dplyr)
library(mlgts)
library(tidyverse)
library(scene)
library(openxlsx)
library(wesanderson)
library(extrafont)
library(ggforce)

source("Q:/VIMC_malaria/VIMC_functions.R", echo=TRUE)

# directories ------------------------------------------------------------------
drat::addRepo("malaria", "file:///f:/drat")
code_dir<- 'Q:/VIMC_malaria/' #  directory where code is stored
malaria_dir<- 'Q:/VIMC_files'      #  project directory where files are stored
setwd('Q:/')

intvn<-read.csv(paste0(malaria_dir, '/output/central_burden_estimates/central_burden_vaccine.csv'))
bl<- read.csv(paste0(malaria_dir, '/output/central_burden_estimates/central_burden_baseline.csv'))

# plot outputs over time  ------------------------------------------------------
intvn<- intvn |> mutate( scenario = 'intervention')
bl<- bl |> mutate(scenario = 'baseline')

output<- rbind(intvn, bl, fill= T)
#font_import()
loadfonts(device = 'win')

pdf(paste0(malaria_dir, '/plots/diagnostics.pdf'))


# load in final outputs from rfp 
output_rfp_bl<- read.xlsx('Q:/VIMC_malaria_rfp/final/central_burden_baseline.xlsx')

output_rfp_bl<- output_rfp_bl |> mutate(scenario = 'baseline')

output_rfp_int<- read.xlsx('Q:/VIMC_malaria_rfp/final/central_burden_vaccine.xlsx')
output_rfp_int<- output_rfp_int |> mutate( scenario = 'intervention')


output_rfp<- rbind(output_rfp_bl, output_rfp_int)


# clinical cases  --------------------------------------------------------------
ggplot(data= output, mapping = aes(x= year, y= cases, color= scenario, fill= scenario))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age, scales= 'free') +
  labs(x= 'Time (in years)', y= 'Clinical cases', title= paste0('Clinical cases over time: ', unique(output$country)),
       color= 'Scenario', fill= 'Scenario') +
  theme_minimal()+
  theme(text= element_text(family= 'Arial')) +
  scale_color_manual(values= wes_palette('Royal2', n= 2)) +
  scale_fill_manual(values= wes_palette('Royal2', n= 2)) 


ggplot(data= output_rfp, mapping = aes(x= year, y= cases/cohort_size, color= scenario, fill= scenario))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap_paginate(~age, 
                      scales = 'free',
                      nrow= 5,
                      ncol= 5) +
  labs(x= 'Time (in years)', y= 'Clinical cases', title= paste0('Clinical cases over time: ', unique(output$country)),
       color= 'Scenario', fill= 'Scenario') +
  theme_minimal()+
  theme(text= element_text(family= 'Arial')) +
  scale_color_manual(values= wes_palette('Royal2', n= 2)) +
  scale_fill_manual(values= wes_palette('Royal2', n= 2)) 

# deaths  ----------------------------------------------------------------------
ggplot(data= output, mapping = aes(x= year, y= deaths, color= scenario, fill= scenario))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~ age) +
  labs(x= 'Time (in years)', y= 'Deaths', title= paste0('Deaths over time: ', unique(output$country)),
       color= 'Scenario', fill= 'Scenario') +
  theme_minimal()+
  theme(text= element_text(family= 'Arial')) +
  scale_color_manual(values= wes_palette('Royal2', n= 2)) +
  scale_fill_manual(values= wes_palette('Royal2', n= 2)) 


# DALYs ------------------------------------------------------------------------
ggplot(data= output, mapping = aes(x= year, y= dalys, color= scenario, fill= scenario))+
  geom_smooth(alpha= 0.2)  +
  facet_wrap(~age) +
  labs(x= 'Time (in years)', y= 'DALYs', title= paste0('DALYs over time: ',unique(output$country)),
       color= 'Scenario', fill= 'Scenario') +
  theme_minimal()+
  theme(text= element_text(family= 'Arial')) +
  scale_color_manual(values= wes_palette('Royal2', n= 2)) +
  scale_fill_manual(values= wes_palette('Royal2', n= 2)) 

dev.off()
