##  title   compare_model_results
##  author  Lydia Haile
##  purpose compare model outputs from RFP and test case

# setup  -----------------------------------------------------------------------
# packages  --------------------------------------------------------------------
library(data.table)
library(openxlsx)
library(dplyr)
library(wesanderson)
library(ggplot2)
library(ggforce)
library(readr)

# directories ------------------------------------------------------------------
dir<- 'Q:/VIMC_malaria_rfp/' #directory issue
save_dir<- 'Q:/VIMC_files/central_estimates/'
source('Q:/VIMC_malaria/VIMC_functions.R')

# load objects 
year <- 365
month <- 30

# make a list of filepaths and tags for plotting
inputs <- data.table(
  'identifier' = c('RFP baseline',
                   'RFP vaccine',
                   'PFPR 17% baseline',
                   'PFPR 17% vaccine',
                   'PFPR 10% baseline',
                   'PFPR 10% vaccine',
                   'PFPR 50% baseline',
                   'PFPR 50% vaccine'),
  'filepath' = c(
    paste0(dir, 'output/central_burden_baseline_0_download.csv'),
    paste0(dir, 'output/central_burden_vaccine_0_download.csv'),
    paste0(save_dir, 'baseline/raw_model_output/raw_model_output_rfp_large_pop_burnin.rds'),
    paste0(save_dir, 'intervention/raw_model_output/raw_model_output_rfp_large_pop_burnin.rds'),
    paste0(save_dir, 'baseline/raw_model_output/raw_model_output/raw_model_output_rfp_PFPR_10.rds'),
    paste0(save_dir, 'intervention/raw_model_output/raw_model_output_rfp_PFPR_10.rds'),
    paste0(save_dir, 'baseline/raw_model_output/raw_model_output_rfp_PFPR_50.rds'),
    paste0(save_dir, 'intervention/raw_model_output/raw_model_output_rfp_PFPR_50.rds'))
)

# format inputs ------------------------------------------------------
#' @title format_for_plots
#' @description
#' Format data for plotting (different model runs)
#'
#' @param input1 identifier for input one
#' @param input2 identifier for input two
#' @param pg     page of the plot to print

#' @export 

  format_for_plots<- function (input1, input2){
    
    filepath1<- inputs[identifier== input1, filepath]
    filepath2<- inputs[identifier== input2, filepath]
    
    
    if (input1 %like% 'RFP') {
      dt1 <- read.csv(filepath1)
      dt1 <- dt1 |>
        rename(clinical = clin_inc)
      
    } else{
      dt1 <- read_rds(filepath1)
    }
    
    if (input2 %like% 'RFP') {
      dt2 <- read.csv(filepath2)
      dt2 <- dt2 |>
        rename(clinical = clin_inc)
      
    } else {
      dt2 <- read_rds(filepath2)
    }
    
    # format inputs for plotting  ----------------------------------------------
    dt1 <- dt1 |>
      mutate(identifier = input1)
    dt2 <- dt2 |>
      mutate(identifier = input2)
    
    dt<- data.table(rbind(dt1, dt2, fill= T))
    
    if(input1 %like% 'PFPR' | input2 %like% 'PFPR'){
      
      dt<- drop_burnin(dt, 15 * 365)
      dt <- dt(
        dt,
        time_divisor = 365,
        baseline_t = 0,
        age_divisor = 1,
        scaler = 0.215,
        treatment_scaler = 0.5,
        baseline_treatment = 0
      )
      
      dt<- dt |>
        rename(year = t) |>
        mutate(age = as.numeric(age_lower)) |>
        mutate(age = round(age/365)) |>
        mutate(year= year + 1999) |>
        select(identifier, year, age, clinical)
      
    }else{
      dt<- dt |>
        mutate(age = as.numeric(age)) |>
        select(identifier, year, age, clinical)
    }
   
    dt<- dt[identifier!= TRUE]
    return(dt) 
  }


# test function to see if it works
input1<- inputs[1]$identifier
input2<- inputs[2]$identifier


#' @title Plot_comparison
#' @description
#' Generate comparison plots for  different model runs
#'
#' @param pg     page of the plot to print
#' @export 

  plot_comparison<- function(dt, pg){
    
    p<- ggplot(data= dt[year < 2051 & age < 70], mapping = aes(x= year, y= clinical, color= identifier, fill= identifier))+
      geom_point(alpha= 0.5)  +
      facet_wrap_paginate(~round(age), scales= 'free', nrow= 5, ncol= 5, page= pg) +
      labs(x= 'Time (in years)', y= 'Clinical incidence rate', 
           title= paste0('Clinical incidence rate over time: (page ', pg, ')'),
           color= 'Model Run', fill= 'Model Run') +
      theme_minimal()+
      scale_color_manual(values= wes_palette('Royal2', n= 2)) +
      scale_x_continuous(guide = guide_axis(n.dodge = 2))+
      scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
    
    
    print(p)
  }
 
  
  

pdf('Q:/VIMC_files/plots/rfp_comparison_plots_baseline_larger_pop.pdf', width= 12, height= 10)

dt<- format_for_plots(input1=  inputs[1]$identifier,
                      input2=  inputs[2]$identifier) 

plot_comparison(dt, pg= 1)
plot_comparison(dt, pg= 2)
plot_comparison(dt, pg= 3)

dev.off()

# make cases averted plots  -----------------------------------------------------------------------------------
inputs <- data.table(
  'identifier' = c('RFP baseline',
                   'RFP vaccine',
                   'PFPR 17% baseline',
                   'PFPR 17% vaccine',
                   'PFPR 10% baseline',
                   'PFPR 10% vaccine',
                   'PFPR 50% baseline',
                   'PFPR 50% vaccine'),
  'filepath' = c(
    paste0(dir, 'final/submitted_estimates/central_burden_baseline.csv'),
    paste0(dir, 'final/submitted_estimates/central_burden_vaccine.csv'),
    paste0(save_dir, 'baseline/raw_model_output_rfp_large_pop_burnin.rds'),
    paste0(save_dir, 'intervention/raw_model_output/raw_model_output_rfp_large_pop_burnin.rds'),
    paste0(save_dir, 'baseline/raw_model_output/raw_model_output_rfp_PFPR_10.rds'),
    paste0(save_dir, 'intervention/raw_model_output/raw_model_output_rfp_PFPR_10.rds'),
    paste0(save_dir, 'baseline/raw_model_output/raw_model_output_rfp_PFPR_50.rds'),
    paste0(save_dir, 'intervention/raw_model_output/raw_model_output_rfp_PFPR_50.rds')))


################################################################################################################

#' @title Cases averted plots
#' @description
#' output: Generate cases averted plots
#'
#' @param baseline             identifier for baseline data table
#' @param intervention         identifier for intervention data table
#' @param tag                  tag for the plot title
#' @export 
#' 

  cases_averted_plot<- function(baseline, intervention, tag){
    
    bl_filepath<- inputs[identifier== baseline, filepath]
    int_filepath<- inputs[identifier== intervention, filepath]
    
    
    if (baseline %like% 'RFP' | intervention %like% 'RFP') {
      # for RFP use cases submitted to VIMC
      
      bl <- read.csv(bl_filepath)
      int <- read.csv(int_filepath)
      
      bl <- bl |>
        rename(cases_bl = cases)
      int<- int |>
        rename(cases_int = cases)
      
      dt <- merge(bl, int, by = c('age', 'year'))
      
      dt <- dt |>
        mutate(cases_averted = cases_bl - cases_int)
      
      dt<- data.table(dt)
      dt<- dt[, cases_averted:= sum(cases_averted), by= c('age')]
      dt<- unique(dt, by= c('age'))
      
    } else{
      
      # for test case merge on population to calculate cases
      pop<- read.xlsx('Q:/VIMC_malaria_rfp/final/central_burden_baseline.xlsx')
      pop<- data.table(pop)
      pop<- pop[, c('age', 'year', 'cohort_size')]
      
      bl <- read_rds(bl_filepath)
      int <- read_rds(int_filepath)
      
      bl<- drop_burnin(bl, 15 * 365)
      int<- drop_burnin(int, 15 * 365)
      
      bl <- get_rates(
        bl,
        time_divisor = 365,
        baseline_t = 0,
        age_divisor = 1,
        scaler = 0.215,
        treatment_scaler = 0.5,
        baseline_treatment = 0
      )
      
      int <- get_rates(
        int,
        time_divisor = 365,
        baseline_t = 0,
        age_divisor = 1,
        scaler = 0.215,
        treatment_scaler = 0.5,
        baseline_treatment = 0
      )
      
      bl<- bl |>
        rename(year = t) |>
        mutate(identifier = baseline) |>
        rename(clinical_bl = clinical) |>
        mutate(age = as.numeric(age_lower)) |>
        mutate(age = round(age/365)) |>
        mutate(year= year + 1999) |>
        select(identifier, year, age, clinical_bl)
      

      int<- int |>
        rename(year = t) |>
        mutate(identifier = intervention) |>
        rename(clinical_int = clinical) |>
        mutate(age = as.numeric(age_lower)) |>
        mutate(age = round(age/365)) |>
        mutate(year= year + 1999) |>
        select(identifier, year, age, clinical_int)

      
      dt <-
        merge(bl, int, by = c('age', 'year'))
      
      dt <-
        merge(dt, pop, by = c('age', 'year'))
      
      dt <- dt |>
        mutate(cases_bl = cohort_size * clinical_bl,
               cases_int = cohort_size * clinical_int)
      
      dt <- dt |>
        mutate(cases_averted = cases_bl - cases_int)
      
      dt<- data.table(dt)
      dt<- dt[year > 2039]
      dt<- dt[, cases_averted:= sum(cases_averted), by= c('age')]
      dt<- unique(dt, by= c('age'))
    }
  
  # generate plot
  p<- ggplot(dt[age < 16], mapping= aes(x= age, y= cases_averted/1000))+
    geom_bar(stat= 'identity', fill= '#2556B896')+
    labs(y= 'Cases averted (in thousands)',
         x= 'Age (in years)',
         title= paste0('Cases averted by age, 2040-2050: ', tag))
  
  
  return(p)
  }

# test function on inputs of different PFPR
  baseline<- inputs[5]$identifier
  intervention<- inputs[6]$identifier
  
  pdf('Q:/VIMC_files/plots/cases_averted_plots_update.pdf', width= 12, height= 10)
  cases_averted_plot(baseline= inputs[7]$identifier,
                     intervention= inputs[8]$identifier,
                     tag= 'PFPR 50%')

  cases_averted_plot(baseline= inputs[5]$identifier,
                     intervention= inputs[6]$identifier,
                     tag= 'PFPR 10%')
  
  dev.off()
  