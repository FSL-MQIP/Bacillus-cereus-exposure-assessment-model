setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model/Re-fit")

library(biogrowth)
library(tibble)
library(EnvStats) # to load rtri function 
library(truncnorm) # to load rtruncnorm function
library(jmuOutlier) # to load rlaplace function

# trace(biogrowth:::secondary_model_data, edit = T)
# change the function to the following codes (add reducedRatkowsky as a new sec model): 
secondary_model_data <- function (model_name = NULL) {
  model_data <- list(CPM = list(identifier = "CPM", 
                                name = "Cardinal Parameter Model", 
                                pars = c("xmin", "xopt", "xmax", "n"), 
                                model = CPM_model, 
                                ref = paste("Rosso, L., Lobry, J. R., Bajard, S., and Flandrois, J. P. (1995).", 
                                      "Convenient Model To Describe the Combined Effects of Temperature and pH on", 
                                      "Microbial Growth. Applied and Environmental Microbiology, 61(2), 610-616.")), 
                     
                     Zwietering = list(identifier = "Zwietering", 
                                       name = "Zwietering gamma function", 
                                       pars = c("xmin", "xopt", "n"), 
                                       model = zwietering_gamma, 
                                       ref = paste("Zwietering, Marcel H., Wijtzes, T., De Wit, J. C., and Riet,", 
                                             "K. V. (1992). A Decision Support System for Prediction of the Microbial", 
                                             "Spoilage in Foods. Journal of Food Protection, 55(12), 973-979.", 
                                             "https://doi.org/10.4315/0362-028X-55.12.973")), 
                     
                     fullRatkowsky = list(identifier = "fullRatkowsky", 
                                          name = "(Adapted) Full Ratkowsky model", 
                                          pars = c("xmin", "xmax", "c"), 
                                          model = full_Ratkowski, 
                                          ref = paste("Ratkowsky, D. A., Lowry, R. K., McMeekin, T. A.,", 
                                                      "Stokes, A. N., and Chandler, R. E. (1983). Model for", 
                                                      "bacterial culture growth rate throughout the entire", 
                                                      "biokinetic temperature range. Journal of Bacteriology,", 
                                                      "154(3), 1222-1226.")),
                     
                     reducedRatkowsky = list(identifier = "reducedRatkowsky", 
                                             name = "Reduced Ratkowsky model", 
                                             pars = c("xmin", "b"),
                                             model = reduced_Ratkowski, 
                                             ref = paste("Ratkowsky, D. A., Lowry, R. K., McMeekin, T. A.,", 
                                            "Stokes, A. N., and Chandler, R. E. (1983). Model for", 
                                            "bacterial culture growth rate throughout the entire", 
                                            "biokinetic temperature range. Journal of Bacteriology,", 
                                            "154(3), 1222-1226.")))
  
                    
  if (is.null(model_name)) {
    return(names(model_data))
  }
  my_model <- model_data[[model_name]]
  if (is.null(my_model)) {
    stop(paste("Unknown model name:", model_name))
  }
  else {
    my_model
  }
}


# trace(biogrowth:::calculate_gammas, edit = T)
# change the function to the following codes (add a line for reduced_Ratkowski):

function (this_t, env_func, sec_models) 
{
  out <- lapply(names(sec_models), function(this_condition) {
    this_x <- env_func[[this_condition]](this_t)
    this_sec <- sec_models[[this_condition]]
    this_gamma <- switch(this_sec$model, 
                         fullRatkowsky = full_Ratkowski(this_x, this_sec$xmin, this_sec$xmax, this_sec$c), 
                         CPM = CPM_model(this_x, this_sec$xmin, this_sec$xopt, this_sec$xmax, this_sec$n), 
                         Zwietering = zwietering_gamma(this_x, this_sec$xmin, this_sec$xopt, this_sec$n), 
                         reducedRatkowsky = reduced_Ratkowski(this_x, this_sec$xmin, this_sec$b), 
                         stop(paste("Model",this_sec$model, "not known.")))
    this_gamma
  })
  out <- unlist(out)
  names(out) <- names(sec_models)
  out
}

# select different xopt by different group names 
xopt_func <- function(group_name){
  if(group_name == "II")
  {return(36.31)}
  else if(group_name == "III")
  {return(39.27)}
  else if(group_name == "IV")
  {return(38.735)}
  else if(group_name == "V")
  {return(37.375)}
  else (group_name == "VII")
  {return(42.35)}
}

## Set up dataframe for modeling 100 units of HTST milk 
n_sim = 100
data = data.frame(unit_id = rep(seq(1,n_sim)))

## Set seed
set.seed(1)

# Stage 1: facility storage 
## (a)  Sample the temperature distribution
data$T_F <- rep(runif(n_sim,min=3.5,max=4.5)) #uniform distribution
## (b) Sample the storage time (in days) distribution
data$t_F <- rep(runif(n_sim,min=1,max=2)) #uniform distribution

# Stage 2: transport from facility to retail store
## (a)  Sample the temperature distribution
data$T_T <- rep(rtri(n_sim,min=1.7,max=10.0,mode=4.4)) #triangular distribution
## (b) Sample the transportation time (in days) distribution
data$t_T <- rep(rtri(n_sim,min=1,max=10,mode=5))

# Stage 3: storage/display at retail store
## (a)  Sample the temperature distribution
data$T_S <- rep(rtruncnorm(n_sim,a=-1.4,b=5.4,mean=2.3,sd=1.8)) #truncated normal distribution
## (b) Sample the storage time (in days) distribution
data$t_S <- rep(rtruncnorm(n_sim,a=0.042,b=10.0, mean=1.821,sd=3.3)) #truncated normal distribution

## Stage 4: transportation from retail store to home
## (a)  Sample the temperature distribution
data$T_T2 <- rep(rtruncnorm(n_sim,a=0,b=10,mean=8.5,sd=1.0)) #truncated normal distribution
## (b) Sample the transportation time (in days) distribution 
data$t_T2 <- rep(rtruncnorm(n_sim,a=0.01,b=0.24, mean=0.04,sd=0.02)) #truncated normal distribution

## Stage 5: home storage 
## (a)  Sample the temperature distribution
temps <- rep(NA, n_sim)
for (i in 1:n_sim){
  number <- rlaplace(1,m=4.06,s=2.31)
  while (number > 15 | number < -1) {
    number <- rlaplace(1,m=4.06,s=2.31) #truncated laplace distribution 
  }
  temps[i] <- number
}
data$T_H <- temps
## (b) Define t_H as 14, 35 days for all units
data$t_H <- rep(35, each = n_sim)

## Model temperature profiles of 100 units HTST milk 
env_cond_time <- matrix(c(rep(0,100),
                          data$t_F, 
                          data$t_F+0.001,
                          data$t_F + data$t_T,
                          data$t_F + data$t_T+0.001,
                          data$t_F + data$t_T + data$t_S,
                          data$t_F + data$t_T + data$t_S+0.001,
                          data$t_F + data$t_T + data$t_S + data$t_T2,
                          data$t_F + data$t_T + data$t_S + data$t_T2+0.001,
                          data$t_F + data$t_T + data$t_S + data$t_T2 + data$t_H), ncol = 10)

env_cond_temp <- matrix(c(data$T_F, 
                          data$T_F,
                          data$T_T,
                          data$T_T,
                          data$T_S,
                          data$T_S,
                          data$T_T2,
                          data$T_T2,
                          data$T_H,
                          data$T_H), ncol = 10)

# define the new secondary model 
reduced_Ratkowski = function(x, xmin, b){
  xopt = xopt_func("IV")                   # change for different B cereus groups to give different xopt                          
  mu_opt = b * (xopt - xmin)
  gamma = b * (x - xmin)
  gamma <- gamma/mu_opt
  gamma <- gamma^2
  gamma[x < xmin] <- 0
  return(gamma)
}

#################################### run the codes below to implement the model####################################
# Import growth parameters
simulation_input = read.csv("InputFiles/simulation_input.csv")
colnames(simulation_input) = c("isolate","Q0","Nmax","b","Tmin","group")
simulation_input$Q0 = as.numeric (simulation_input$Q0)
simulation_input$Nmax = as.numeric(simulation_input$Nmax)
simulation_input$Tmin = as.numeric(simulation_input$Tmin)
simulation_input$b = as.numeric(simulation_input$b)

# Subset data for each isolate 
Iso649 = subset(simulation_input,isolate == 649)

# Prepare model input
xmin = Iso649$Tmin
b = Iso649$b 
mu_opt = (b*(xopt_func("IV")-xmin))^2 # assume (Topt,sqrt(mu_opt)) is on the linear region 
Q0 = Iso649$Q0

# Run the models 
my_primary <- list(mu_opt = mu_opt, 
                   Nmax = Iso649$Nmax,     
                   N0 = 1e2,          # depending on product testing result, assumed in this case 
                   Q0 = Q0)

sec_temperature <- list(model = "reducedRatkowsky",  
                        xmin = xmin, 
                        b = b)    

my_secondary <- list(temperature = sec_temperature)

for (i in 1:n_sim){
  growth = predict_dynamic_growth(times = env_cond_time[i,],
                                        env_conditions = tibble(time = env_cond_time[i,],
                                                                temperature = env_cond_temp[i,]),
                                        my_primary,
                                        my_secondary)
  sim = growth$simulation
  data$conc[i] = tail(sim$logN, 1)
}

# Generate output 
data_649_d35<-data
Iso649_d35_5log = sum(data_649_d35$conc>5)/100

simulation_result_d14 <- rbind(data_193_d14,data_194_d14,data_402_d14,data_407_d14,
                    data_413_d14,data_433_d14,data_457_d14,data_474_d14, 
                    data_495_d14,data_518_d14,data_536_d14,data_564_d14,
                    data_570_d14,data_638_d14,data_649_d14)

simulation_result_d35 <- rbind(data_193_d35,data_194_d35,data_402_d35,data_407_d35,
                               data_413_d35,data_433_d35,data_457_d35,data_474_d35, 
                               data_495_d35,data_518_d35,data_536_d35,data_564_d35,
                               data_570_d35,data_638_d35,data_649_d35)

prediction_result_d14 <- c(0.03,0.03,0.05,0.03,0.03,0.02,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.04,0.03)
prediction_result_d35 <- c(0.04,0.05,0.15,0.03,0.13,0.04,0.06,0.05,0.04,0.05,0.04,0.1,0.1,0.13,0.04)
prediction_result <- cbind(prediction_result_d14,prediction_result_d35)
rownames(prediction_result) <- c("193","194","402","407","413","433","457","474",
                                 "495","518","536","564","570","638","649")

mean_d14 = mean(prediction_result_d14)
sd_d14 = sd(prediction_result_d14)

mean_d35 = mean(prediction_result_d35)
sd_d35 = sd(prediction_result_d35)

mean = c(0.03,0.07)
sd = c(0.0064,0.0404)
prediction_result_refit = rbind(prediction_result,mean,sd)
write.csv(prediction_result_refit,"OutputFiles/prediction_result_refit.csv")
