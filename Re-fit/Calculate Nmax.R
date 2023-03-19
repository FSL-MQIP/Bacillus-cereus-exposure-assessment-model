# Set work directory 
setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model/Re-fit")

# Import data set
Nmax_table <- read.csv("InputFiles/Nmax.csv")
colnames(Nmax_table) <- c("isolate", "rep", "LOG10Nmax","T")

# calculate average Nmax by isolate
avg_Nmax_by_isolate <- aggregate(LOG10Nmax ~ isolate, data = Nmax_table, mean)
colnames(avg_Nmax_by_isolate) <- c("isolate", "average_LOG10Nmax")
avg_Nmax_by_isolate$average_Nmax = 10^avg_Nmax_by_isolate$average_LOG10Nmax

# export the data to a CSV file
write.csv(avg_Nmax_by_isolate, "OutputFiles/Nmax_new.csv", row.names = FALSE)


