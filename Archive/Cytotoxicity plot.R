library(tidyr)

data = read.csv("FINAL B_cereus_mastersheet.csv")
colnames(data)[1] <- "Isolate"
data <- data[, c("Isolate", "Btyper3.Closest.Type.Strain", "Adjusted.panC.Group..predicted.species.","Average.Cell.Viability")]

data <- separate(data, col = Btyper3.Closest.Type.Strain, into = c("species", "ANI"), sep = "\\(")
data <- separate(data, col = Adjusted.panC.Group..predicted.species., into = c("clade", "Predicted.species"), sep = "\\(")
data = data[,c("Isolate","species","clade","Average.Cell.Viability")]
colnames(data) [4] <- "Normalized_Cytotoxicity"

# Match clade and species, display histogram of all isolates, matched clade and species 
min_value <- min(data$Normalized_Cytotoxicity)
max_value <- max(data$Normalized_Cytotoxicity)
breaks <- seq(floor(min_value), ceiling(max_value) + 0.1, by = 0.1)
hist(data$Normalized_Cytotoxicity, 
     breaks = breaks, 
     xlab = "Normalized_Cytotoxicity", 
     ylab = "Frequency",
     main = "Histogram of Normalized Cytotoxicity")

data_groupII <- subset(data, clade == "Group_II")
hist(data_groupII$Normalized_Cytotoxicity, breaks = breaks, add = TRUE, col = "yellow", border = "black")

data_tropicus <- subset(data, species == "tropicus")
hist(data_tropicus$Normalized_Cytotoxicity, breaks = breaks, add = TRUE, col = "blue", border = "black")

abline(v = 0.3, col = "red", lty = "dashed", lwd = 2)

legend("topright", 
       legend = c("All isolates", "clade II", "tropicus"), 
       fill = c("gray", "yellow","blue"), border = "black")
