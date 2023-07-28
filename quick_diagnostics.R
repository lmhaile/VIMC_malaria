#quick diagnostics


filepath<- files[[1]]
dt<- readRDS(filepath)
raw<- copy(dt)

# raw model output
ggplot(data= raw, mapping= aes(x= timestep, y= n_0_364))+
  geom_point() +
  #facet_wrap(~age) +
  labs(y= 'Population',
       x= 'year',
       title= paste0('Population over time by age: ', site_name, ' ', ur )) +
  theme_minimal()


# quick diagnostic
# test<- dt |>
#   filter(t == 2005)
# 
# sum(test$prop_n)
ggplot(data= pop , mapping= aes(x= t, y= pop))+
  geom_point() +
  labs(y= 'Population',
       x= 'year',
       title= paste0('Population over time: ', site_name, ' ', ur)) +
  theme_minimal()



ggplot(data= dt, mapping= aes(x= timestep/365, y= population))+
  geom_point() +
  facet_wrap(~age_years_lower) +
  labs(y= 'Population',
       x= 'year',
       title= paste0('Population over time by age: ', site_name, ' ', ur )) +
  theme_minimal()
