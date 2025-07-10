# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Sample data (replace this with your actual data reading process)
data <- data.frame(
  Gene = c("CYB5A", "CYB5A", "CYB5A", "CYB5A", "CYB5A", "CYB5A", "CYP17A1", "CYP17A1", "SULT2A1", "SULT2A1"),
  Chr = c(18, 18, 18, 18, 18, 18, 10, 10, 19, 19),
  Position = c(71936427, 71943304, 71947422, 71948446, 71955452, 71958555, 104597132, 104597159, 48388658, 48389363),
  AFR = c(0, 0, 0, 0, 0, 0, 0, 0.001, 0.001, 0.001),
  EAS = c(0.070, 0.061, 0.059, 0.059, 0.058, 0.057, 0.006, 0.000, 0.003, 0.004),
  EUR = c(0.112, 0.112, 0.113, 0.113, 0.112, 0.112, 0.010, 0.000, 0.133, 0.135),
  AMR = c(0.217, 0.223, 0.219, 0.212, 0.209, 0.210, 0.000, 0.000, 0.059, 0.059),
  SAS = c(0.105, 0.104, 0.104, 0.104, 0.104, 0.104, 0.000, 0.000, 0.087, 0.088)
) #each data point in the superpop parentheses are the superpop average for the specific position

data_long <- data %>%
  pivot_longer(cols = c(AFR, EAS, EUR, AMR, SAS), 
               names_to = "Superpopulation", 
               values_to = "Genetic_Distance") %>%
  # Create a new column for combined chromosome, position, and gene labels
  mutate(Label = paste0(Gene, " (", Chr, ":", Position, ")"))

ggplot(data_long, aes(x = Superpopulation, y = reorder(Label, Position), fill = Genetic_Distance)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "blue", high = "red", name = "Alternate Allele") + 
  labs(title = "Heatmap of Neanderthal Alternate Allele Across Superpopulations", 
       x = "Superpopulation", 
       y = "Genomic Position") +  # Removed the extra "+" here
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 8))

