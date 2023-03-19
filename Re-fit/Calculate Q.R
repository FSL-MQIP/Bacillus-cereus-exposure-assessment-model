# Set work directory 
setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model/Re-fit")

# Import data set
h0_table <- read.csv("InputFiles/h0 table.csv")
colnames(h0_table) <- c("isolate","rep","T","lag(h)","lag(day)","mumax (ln CFU/mL per h)","mumax (log10 CFU/mL per day)","h0")  

# Calculate Q0 for each isolate
h0_table$Q0 <- 1/(exp(h0_table$h0)-1)

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
