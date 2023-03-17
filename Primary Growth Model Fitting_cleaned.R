# load packages
library(nlsMicrobio)
library(minpack.lm)

# load primary growth models
source("UtilityFunctions_baranyi.R")

# fit 22 dC growth data 
# import growth data set 
# 22dC growth data 
data1 <-read.csv("InputFiles/data_22_new.csv")

# create a list of data frames, one for each isolate and rep
data_list <- split(data1, list(data1$Isolate, data1$Rep))

# 22dC starting values
starting_values_22dC <-read.csv("InputFiles/Starting values_22dC.csv")

# subset data and fit Baranyi model
# initialize a list to store the fits
fit_list <- list()

# loop over each sample and its repetitions
for (i in 1:length(data_list)) {
  # subset the data for the current sample and rep
  sub_data <- data_list[[i]]
  
  # extract the sample and rep information
  isolate <- unique(sub_data$Isolate)
  rep <- unique(sub_data$Rep)
  
  # set the starting values for the current isolate and rep
  start_values <- c(LOG10N0 = starting_values_22dC$LOG10N0[i],
                    lag = starting_values_22dC$lag[i],
                    mumax = starting_values_22dC$mumax[i],
                    LOG10Nmax = starting_values_22dC$LOG10Nmax[i])
  
  # fit the Baranyi model to the subset of data
  fit <- nlsLM(LOG10N ~ baranyi_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=sub_data,
               start = start_values, 
               lower = c(0,0,0,0))
  
  # add the fit and its summary to the fit list
  fit_list[[i]] <- data.frame(isolate=isolate, rep=rep, coef=coefficients(fit))
  fit_list[[i]]$rsq <- summary(fit)$r.squared
}

# combine data and generate output 
# combine the fits into a single data frame
fits <- do.call(rbind, fit_list)

# generate output
write.csv(fits,"gp_22dC_new.csv")

