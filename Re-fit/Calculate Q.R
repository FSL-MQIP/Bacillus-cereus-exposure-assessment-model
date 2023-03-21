# Set work directory 
setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model/Re-fit")

# Import data set
data_22dC <- read.csv("OutputFiles/gp_22dC_new.csv")
data_22dC$T <- rep(22,34)
data_10dC <- read.csv("OutputFiles/gp_10dC_new.csv")
data_10dC$T <- rep(10,30)

# Generate h0 table
h0_table <- rbind(data_22dC,data_10dC)
h0_table <- h0_table[,c("isolate","rep","T","mumax",'lag')]
h0_table$mumax2 <- h0_table$mumax*24
h0_table$lag2 <- h0_table$lag/24
colnames(h0_table) <- c("isolate","rep","T","mumax_ln_h","lag_h","mumax_ln_day","lag_day")  
h0_table$h0 <- h0_table$mumax_ln_day*h0_table$lag_day
h0_table <- h0_table[order(h0_table$isolate), ]
h0_table <- h0_table[h0_table$lag_h != 0,]

# Save file
write.csv(h0_table,"OutputFiles/h0 table.csv")

# Calculate Q0 for each isolate
h0_table$Q0 <- 1/(exp(h0_table$h0)-1)    # mumax is in natural (i.e. exp(1)) scale to calculate h0 in this formula

# Group the data by isolate and calculate the mean Q0 for each group
grouped_Q0 <- aggregate(h0_table$Q0, by = list(h0_table$isolate), FUN = mean)
colnames(grouped_Q0) <- c("isolate", "Q0")

# Calculate the mean h0 for each isolate
avg_h0 <- aggregate(h0 ~ isolate, data = h0_table, FUN = mean)

# Calculate the standard deviation of h0 for each isolate
sd_h0 <- aggregate(h0 ~ isolate, data = h0_table, FUN = sd)

# Merge the two data frames
grouped_h0 <- merge(avg_h0, sd_h0, by = "isolate")

# Rename the columns
colnames(grouped_h0) <- c("isolate", "mean_h0", "sd_h0")

# Merge the Q0 and h0 data frames
result <- merge(grouped_Q0, grouped_h0, by = "isolate")

# Export the results to a CSV file
write.csv(result, file = "OutputFiles/Q0_h0_summary.csv", row.names = FALSE)
