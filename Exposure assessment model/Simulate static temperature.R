setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model")

# load packages 
library(tibble)
library(dplyr)
library(tidyr)
library(formula.tools)    # to load 'rhs'
library(purrr)            # to load 'map'
library(deSolve)          # to load 'ode'
library(rlang)

# load utility functions
source("Exposure assessment model/UtilityFunctions_dynamic_growth.R")

# predict dynamic growth
# import data set
data_Q0 = read.csv("Exposure assessment model/InputFiles/Q0_h0_summary.csv")
data_Nmax = read.csv("Exposure assessment model/InputFiles/Nmax_new.csv")
data_sec_model = read.csv("Exposure assessment model/InputFiles/sec_model_new.csv")
Clade = c("I","II","VII","IV","IV","IV","IV","II","III","IV","II","VII","II","V","V","IV")

# generate simulation input
simulation_input <- data.frame(isolate = data_Q0$isolate, Q0 = data_Q0$Q0, 
                               Nmax = data_Nmax$average_Nmax, b = data_sec_model$Slope, 
                               Tmin = data_sec_model$Tmin,Clade)
simulation_input$Topt = sapply(simulation_input$Clade, xopt_func)
simulation_input$mu_opt = (simulation_input$b*(simulation_input$Topt-simulation_input$Tmin))^2 

# create a tibble with environmental conditions over time (Tmin)
env_cond_time <- c(0:90)
env_cond_temp <- data_sec_model$Tmin

# run simulation (Tmin)
count_d90_Tmin <- vector("list",16)

for(i in 1:16) {
  my_primary <- list(mu_opt = simulation_input$mu_opt[i], 
                     Nmax = simulation_input$Nmax[i], 
                     N0 = 1e2, 
                     Q0 = simulation_input$Q0[i])
  
  sec_temperature <- list(model = "reducedRatkowsky", 
                          xmin = simulation_input$Tmin[i], 
                          b = simulation_input$b[i], 
                          clade = simulation_input$Clade[i])

  my_secondary <- list(temperature = sec_temperature)

  growth <- predict_dynamic_growth(times = env_cond_time,
                                   env_conditions = tibble(time = env_cond_time,
                                                           temperature = rep(env_cond_temp[i],91)),
                                                           my_primary,
                                                           my_secondary)
  sim <- growth$simulation
  count_d90_Tmin[[i]] = tail(sim$logN, 1)
}
count_d90_Tmin
count_d90_Tmin <- as.data.frame(count_d90_Tmin)
simulation_Tmin_result <- rbind(simulation_input$isolate,count_d90_Tmin)
write.csv(simulation_Tmin_result,"Exposure assessment model/OutputFiles/simulation_Tmin_result.csv")

# create a tibble with environmental conditions over time (7dC)
env_cond_time <- c(0:90)
env_cond_temp <- 7
simulation_input_638 <- subset(simulation_input,isolate == 638)

# run simulation (7dC)
my_primary <- list(mu_opt = simulation_input_638$mu_opt, 
                   Nmax = simulation_input_638$Nmax, 
                   N0 = 1e2, 
                   Q0 = simulation_input_638$Q0)

sec_temperature <- list(model = "reducedRatkowsky", 
                        xmin = simulation_input_638$Tmin, 
                        b = simulation_input_638$b, 
                        clade = simulation_input_638$Clade)

my_secondary <- list(temperature = sec_temperature)

growth <- predict_dynamic_growth(times = env_cond_time,
                                 env_conditions = tibble(time = env_cond_time,
                                                         temperature = rep(env_cond_temp,91)),
                                 my_primary,
                                 my_secondary)
sim <- growth$simulation
count_d90_638 = tail(sim$logN, 1)

# create a tibble with environmental conditions over time (12dC)
env_cond_time <- c(0:90)
env_cond_temp <- 12
simulation_input_125 <- subset(simulation_input,isolate == 125)

# run simulation (12dC)
my_primary <- list(mu_opt = simulation_input_125$mu_opt, 
                   Nmax = simulation_input_125$Nmax, 
                   N0 = 1e2, 
                   Q0 = simulation_input_125$Q0)

sec_temperature <- list(model = "reducedRatkowsky", 
                        xmin = simulation_input_125$Tmin, 
                        b = simulation_input_125$b, 
                        clade = simulation_input_125$Clade)

my_secondary <- list(temperature = sec_temperature)

growth <- predict_dynamic_growth(times = env_cond_time,
                                 env_conditions = tibble(time = env_cond_time,
                                                         temperature = rep(env_cond_temp,91)),
                                 my_primary,
                                 my_secondary)
sim <- growth$simulation
count_d90_125 = tail(sim$logN, 1)

# create a tibble with environmental conditions over time (9dC)
env_cond_time <- c(0:90)
env_cond_temp <- 9

# subset the data to exclude isolates 125 and 638
simulation_input_sub <- simulation_input %>%
  filter(isolate != 125 & isolate != 638)

# group the data by isolate name
simulation_input_groups <- simulation_input_sub %>% 
  group_by(isolate)

# run simulation (9dC)
count_d90_by_isolate <- simulation_input_groups %>%
  nest() %>%
  mutate(
    count_d90 = map(data, ~{
      my_primary <- list(mu_opt = .$mu_opt, Nmax = .$Nmax, N0 = 1e2, Q0 = .$Q0)
      sec_temperature <- list(model = "reducedRatkowsky", xmin = .$Tmin, b = .$b, clade = .$Clade)
      my_secondary <- list(temperature = sec_temperature)
      growth <- predict_dynamic_growth(
        times = env_cond_time,
        env_conditions = tibble(time = env_cond_time, temperature = rep(env_cond_temp, 91)),
        my_primary,
        my_secondary
      )
      sim <- growth$simulation
      tail(sim$logN, 1)
    })
  ) %>%
  select(isolate, count_d90)

count_d90_by_isolate <- as.data.frame(count_d90_by_isolate)
simulation_result = rbind(count_d90_by_isolate,c(125,count_d90_125),c(638,count_d90_638))
simulation_result = do.call(rbind, simulation_result)
write.csv(simulation_result,"Exposure assessment model/OutputFiles/simulation_result.csv")

