######################################################
## title    Visualize results
## author   Lydia Haile
## purpose  functions used to visualize results
######################################################


#  write up model diagnostics for test site  -----------------------------------
#' Plot model PFPR against prevalence 
#'
#' @param   dt               raw model output
#' @returns plot of PFPR2-10 against prevalence data for site
plot_model_against_prevalence<- function(dt){
  
  require(data.table)
  
  site_data<- foresite::get_site(unique(dt$iso))
  
  site_name<- unique(dt$site_name)
  ur<- unique(dt$urban_rural)
  
  prevalence<- as.data.table(site_data$prevalence)
  
  prevalence<- prevalence[name_1== site_name & urban_rural== ur]
  
  dt$prevalence <- dt$n_detect_730_3649 / dt$n_730_3649
  
  
  # Set the time
  dt$year <- (dt$timestep / 365) + 2000
  
  p<- ggplot()+
    geom_line(data= dt[year < 2021], mapping= aes(x= year, y= prevalence), linewidth= 0.25)+
    geom_point(data= prevalence, mapping= aes(x= year, y= pfpr), size= 2, pch= 19)+
    labs(title= 'Model output against parasite prevalence',
         y= 'PfPR',
         x= 'Year') +
    theme_minimal()+
    theme(text= element_text(family= 'Verdana'))+
    scale_color_manual(values= wes_palette('Royal2', n= 2)) +
    scale_fill_manual(values= wes_palette('Royal2', n= 2)) 
  
  return(p)
}




