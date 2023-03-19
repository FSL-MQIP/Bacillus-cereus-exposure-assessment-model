# Set work directory 
setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model/Re-fit")

# Import data set
data_22dC <- read.csv("OutputFiles/gp_22dC_new.csv")
data_22dC$T <- rep(22,34)
data_10dC <- read.csv("OutputFiles/gp_10dC_new.csv")
data_10dC$T <- rep(10,30)

# Generate Nmax table
Nmax_table <- rbind(data_22dC,data_10dC)
Nmax_table <- Nmax_table[,c("isolate","rep","LOG10Nmax","T")]
colnames(Nmax_table) <- c("isolate", "rep", "LOG10Nmax","T")

# calculate average Nmax by isolate
avg_Nmax_by_isolate <- aggregate(LOG10Nmax ~ isolate, data = Nmax_table, mean)
colnames(avg_Nmax_by_isolate) <- c("isolate", "average_LOG10Nmax")
avg_Nmax_by_isolate$average_Nmax = 10^avg_Nmax_by_isolate$average_LOG10Nmax

# export the data to a CSV file
write.csv(avg_Nmax_by_isolate, "OutputFiles/Nmax_new.csv", row.names = FALSE)


