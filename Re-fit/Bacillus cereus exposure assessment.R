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
                                             pars = c("xmin", "b", "clade"),
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
                         reducedRatkowsky = reduced_Ratkowski(this_x, this_sec$xmin, this_sec$b, this_sec$clade), 
                         stop(paste("Model",this_sec$model, "not known.")))
    this_gamma
  })
  out <- unlist(out)
  names(out) <- names(sec_models)
  out
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

## Define function to select Topt by groups (from literature, create citation)
xopt_func <- function(clade){
  if(clade == "II")
  {return(36.31)}
  else if(clade == "III")
  {return(39.27)}
  else if(clade == "IV")
  {return(38.735)}
  else if(clade == "V")
  {return(37.375)}
  else (clade == "VII")
  {return(42.35)}
}

# Import data set
data_Q0 = read.csv("OutputFiles/Q0_h0_summary.csv")
data_Nmax = read.csv("OutputFiles/Nmax_new.csv")
data_sec_model = read.csv("OutputFiles/sec_model_new.csv")
Clade = c("II","VII","IV","IV","IV","IV","II","III","IV","II","VII","II","V","V","IV")

# Generate simulation input
simulation_input <- data.frame(isolate = data_Q0$isolate[3:17], Q0 = data_Q0$Q0[3:17], Nmax = data_Nmax$average_Nmax[3:17], 
                               b = data_sec_model$b[3:17], Tmin = data_sec_model$Tmin[3:17],Clade)
simulation_input$Topt = sapply(simulation_input$Clade, xopt_func)
simulation_input$mu_opt = (simulation_input$b*(simulation_input$Topt-simulation_input$Tmin))^2 

## Define new secondary model 
reduced_Ratkowski = function(x, xmin, b, clade){
  mu_opt = b * (xopt_func(clade) - xmin)
  gamma = b * (x - xmin)
  gamma <- gamma/mu_opt
  gamma <- gamma^2
  gamma[x < xmin] <- 0
  return(gamma)
}

# Run simulation for each isolate 
# initialize a list to store the final concentrations
final_conc_d35 <- vector(mode = "list", length = nrow(simulation_input))

# loop over each sample
for (i in 1:nrow(simulation_input)) {
  
  # initialize a list to store the final concentrations for this isolate
  final_conc_isolate <- vector(mode = "list", length = n_sim)
  
  # subset data by isolate
  Iso <- simulation_input[i,]
  
  # Prepare model input
  xmin = Iso$Tmin
  b = Iso$b 
  mu_opt = Iso$mu_opt # assume (Topt,sqrt(mu_opt)) is on the linear region 
  Q0 = Iso$Q0
  
  # Run the models 
  my_primary <- list(mu_opt = mu_opt,   # in log10 scale
                     Nmax = Iso$Nmax,     
                     N0 = 1e2,          
                     Q0 = Q0)
  
  sec_temperature <- list(model = "reducedRatkowsky",  
                          xmin = xmin, 
                          b = b,
                          clade = Iso$Clade)    
  
  my_secondary <- list(temperature = sec_temperature)
  
  for (j in 1:n_sim){
    growth = predict_dynamic_growth(times = env_cond_time[j,],
                                    env_conditions = tibble(time = env_cond_time[j,],
                                                            temperature = env_cond_temp[j,]),
                                    my_primary,
                                    my_secondary)
    sim = growth$simulation
    
    # store the final concentration for this simulation
    final_conc_isolate[[j]] <- tail(sim$logN, 1)
    
  }
  
  # store the final concentrations for all simulations for this isolate
  final_conc_d35[[i]] <- final_conc_isolate
}

# Convert the list of 15 elements into a matrix
matrix <- t(sapply(final_conc_d35, unlist))

# Calculate the percentage over 5 for each row (element) in the matrix
percent_over_5_d35 <- rowMeans(matrix > 5) * 100
simulation_input$percent_over_5_d35 = percent_over_5_d35

# Output 
Result = simulation_input
write.csv(Result,"OutputFiles/Result.csv")

