setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model/Exposure assessment model/plot")

library(ggplot2)
library(tidyverse)

gp = read.csv("GP.csv")
Result = read.csv("Result.csv")
data <- read.csv("Tornado.csv")
colnames(data)[1] = "isolate"

# Figure 1. Boxplot that shows the within and between group variability of growth parameters
# Filter out Group III from the dataset
gp_filtered <- gp[gp$Group != "III", ]

# Define temperature levels
temperature_levels <- c("22", "10", "16")

# Set the layout to have 2 rows and 3 columns
par(mfrow = c(2, 3))

lighter_green <- rgb(144, 255, 144, maxColorValue = 255)
lighter_blue <- rgb(173, 216, 255, maxColorValue = 255)

# Loop through temperature levels
for (i in 1:3) {
  # Subset the data for the current temperature level
  current_data <- subset(gp_filtered, Temp == temperature_levels[i])
  
  # Create boxplot for lag on the top row
  boxplot(current_data$lag ~ factor(current_data$Group), col = lighter_green,
          xlab = "Group",
          ylab = paste("lag (h)"),
          main = LETTERS[i])
}

# Loop through temperature levels again for mumax on the bottom row
for (i in 1:3) {
  # Subset the data for the current temperature level
  current_data <- subset(gp_filtered, Temp == temperature_levels[i])
  
  # Create boxplot for mumax on the bottom row
  boxplot(current_data$mumax ~ factor(current_data$Group), col = lighter_blue,
          xlab = "Group",
          ylab = paste("µmax (ln/h)"),
          main = LETTERS[i + 3])
}

# Figure 2. Base model prediction results 
orange <- rgb(255, 165, 0, maxColorValue = 255)
green <- rgb(100, 255, 150, maxColorValue = 255)
blue <- rgb(100, 200, 255, maxColorValue = 255)

ggplot(Result, aes(x = as.factor(isolate), y = percent_over_5, fill = factor(day))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(vars(Group), scales = "free") +
  labs(x = "Isolate name", y = "Percent milk containers > 5 log (%)") +
  scale_fill_manual(values = c(orange, blue, green), 
                    labels = c("Day 14", "Day 21", "Day 35"),
                    name = "Consumer Storage Day") + 
                    ylim(0,15) + 
                    ggtitle("A") +
                    theme_minimal()

ggplot(Result, aes(x = as.factor(isolate), y = percent_over_3, fill = factor(day))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(vars(Group), scales = "free") +
  labs(x = "Isolate name", y = "Percent milk containers > 3 log (%)") +
  scale_fill_manual(values = c(orange, blue, green), 
                    labels = c("Day 14", "Day 21", "Day 35"),
                    name = "Consumer Storage Day") +
                    ylim(0,15) + 
                   ggtitle("B") +
                   theme_minimal() 

# Figure 3. Tornado diagram for sensitivity analysis of Q0 
red <- rgb(178, 34, 34, maxColorValue = 255)
blue <- rgb(100, 200, 255, maxColorValue = 255)

data_long <- tidyr::pivot_longer(data, cols = c("Q0_min", "Q0_max"), names_to = "range_type", values_to = "value")
ggplot(data_long, aes(x = factor(isolate, levels = unique(data$isolate)), y = value, fill = range_type)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_text(data = subset(data_long, range_type == "Q0_min"), 
            aes(label = isolate), 
            hjust = 1.5, vjust = 0.5, size = 3) +
  coord_flip() +
  labs(title = NULL,
       x = "Isolate",
       y = "Percentage point increases and decreasse from the base model prediction results") +
  scale_fill_manual(values = c("Q0_min" = blue, "Q0_max" = red)) +
  theme_minimal() +
  theme(axis.text.y = element_blank()) + 
  guides(fill = guide_legend(title = NULL))
