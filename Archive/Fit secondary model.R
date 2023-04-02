# set work directory
setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model")

# import data set
data <- read.csv("InputFiles/mumax_new.csv")

# change column name 
colnames(data) <- c("Isolate","T","mumax..ln.CFU.mL.per.h.","mumax..log10.CFU.mL.per.day.","sqrt_mumax..sqrt_log10.CFU.mL.per.day.")

# subset data by isolates
unique_isolates <- unique(data[c("Isolate")])

# initialize a list to store the fits  
b <- list()
Tmin <- list()

# loop over each isolate and rep combination to fit secondary model
for (i in 1:nrow(unique_isolates)) {
  isolate <- unique_isolates$Isolate[i]
  isolate_data <- subset(data, Isolate == isolate) # subset the data for the current isolate
  sqrt_mumax <- isolate_data$sqrt_mumax..sqrt_log10.CFU.mL.per.day.
  fit <- lm(sqrt_mumax ~ isolate_data$T)
  b[isolate] <- coef(fit)[2]
  Tmin[isolate] <- -coef(fit)[1]/coef(fit)[2]
}

# create a data frame to store the results
results <- data.frame(
  Isolate = character(),
  Slope = numeric(),
  Tmin = numeric(),
  stringsAsFactors = FALSE
)

# loop over each isolate and add the results to the data frame
for (i in 1:nrow(unique_isolates)) {
  isolate <- unique_isolates$Isolate[i]
  slope <- b[[isolate]]
  tmin <- Tmin[[isolate]]
  results <- rbind(results, data.frame(Isolate = isolate, b = slope, Tmin = tmin))
}

# write the results to a CSV file
write.table(results, file = "OutputFiles/sec_model_new.csv", sep = ",", row.names = FALSE, quote = TRUE)



