#### Load each database
library(data.table)
library(dplyr)
data <- fread('./Training_Set.csv', header = T, stringsAsFactors = F)

## Remove Null Data
# Data
data <- na.omit(data)
data <- data[data$Ikey != "0" & data$Ikey != "#N/A" & data$Ikey != "" & data$Ikey != "NA" & data$Ikey != "NaN", ]
sum(!grepl("-", data$Ikey)) # [1] 0
data <- data[data$UNIPROT_AC != "0" & data$UNIPROT_AC != "#N/A" & data$UNIPROT_AC != "" & data$UNIPROT_AC != "NA" & data$UNIPROT_AC != "NaN", ]
data <- data[!grepl("|", data$UNIPROT_AC, fixed = T), ] # fixed = T: Exact match

#### Load library
library(tidyr)
library(ggplot2)
library(caret)
library(Metrics)

refined_data <- data
### Based on ligand
# Sort train, validation, and test sets by Ikey
refined_data <- refined_data %>% arrange(Ikey)

## Calculate average binding affinity for each Ikey in training data
average_affinity <- refined_data %>%
  group_by(Ikey) %>%
  summarise(Avg_pValue = mean(pValue)) %>%
  arrange(Avg_pValue)
average_affinity <- average_affinity %>% mutate(Idx = row_number())

# Merge index back into the original training data
refined_data <- refined_data %>%
  left_join(average_affinity, by = "Ikey") %>%
  dplyr::select(Ikey, UNIPROT_AC, pValue, Idx)

# Split data into initial training (80%) and test (20%)
set.seed(123) # For reproducibility
initial_split <- createDataPartition(refined_data$pValue, p = 0.8, list = FALSE)
train_data <- refined_data[initial_split, ]
test_data <- refined_data[-initial_split, ]

# Polynomial Regression
poly_model <- lm(pValue ~ poly(Idx, degree = 3), data = train_data)

# Make predictions on the test set
test_predictions_poly <- predict(poly_model, test_data)
pcc_poly_test_cmp <- cor(test_data$pValue, test_predictions_poly)
test_data$Predicted <- test_predictions_poly

# Create the ggplot with test data points
test_data <- test_data %>% filter(pValue > -9)
set.seed(123)

test_data_cmp <- test_data
p1 <- ggplot(test_data_cmp, aes(x = Idx, y = pValue)) +
  geom_point(size = 1) +  # Add the data points as dots
  geom_line(aes(y = Predicted), color = "red") +  # Add the prediction line from the model
  annotate("text", x = mean(test_data_cmp$Idx), y = max(test_data_cmp$pValue, na.rm = TRUE), 
           label = paste0("Cubic regression (PCC = ", round(pcc_poly_test_cmp, 3), ")"), 
           hjust = 0.4, vjust = 1.25, color = "red", size = 5) +  # Add the R-squared value
  labs(x = "Rank compound index", y = "Compound-protein binding affinity") +  # Label the axes
  theme_minimal() + # Use a minimal theme
  theme(axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))


### Based on protein
# Sort train, validation, and test sets by Ikey
refined_data <- data
refined_data <- refined_data %>% arrange(UNIPROT_AC)

## Calculate average binding affinity for each Ikey in training data
average_affinity <- refined_data %>%
  group_by(UNIPROT_AC) %>%
  summarise(Avg_pValue = mean(pValue)) %>%
  arrange(Avg_pValue)
average_affinity <- average_affinity %>% mutate(Idx = row_number())

# Merge index back into the original training data
refined_data <- refined_data %>%
  left_join(average_affinity, by = "UNIPROT_AC") %>%
  dplyr::select(Ikey, UNIPROT_AC, pValue, Idx)

# Split data into initial training (80%) and test (20%)
set.seed(123) # For reproducibility
initial_split <- createDataPartition(refined_data$pValue, p = 0.8, list = FALSE)
train_data <- refined_data[initial_split, ]
test_data <- refined_data[-initial_split, ]

# Polynomial Regression
poly_model <- lm(pValue ~ poly(Idx, degree = 3), data = train_data)

# Make predictions on the test set
test_predictions_poly <- predict(poly_model, test_data)
pcc_poly_test_prn <- cor(test_data$pValue, test_predictions_poly)
test_data$Predicted <- test_predictions_poly

# Create the ggplot with test data points
test_data <- test_data %>% filter(pValue > -10)
test_data_prn <- test_data
p2 <- ggplot(test_data_prn, aes(x = Idx, y = pValue)) +
  geom_point(size = 1) +  # Add the data points as dots
  geom_line(aes(y = Predicted), color = "red") +  # Add the prediction line from the model
  annotate("text", x = max(test_data_prn$Idx), y = min(test_data_prn$pValue, na.rm = TRUE), 
           label = paste0("Cubic regression (PCC = ", round(pcc_poly_test_prn, 3), ")"), # Third-degree Polynomial Regression
           hjust = 1.35, vjust = -2.5, color = "red", size = 5) +  # Add the R-squared value
  labs(x = "Rank protein index", y = "Compound-protein binding affinity") +  # Label the axes
  theme_minimal() + # Use a minimal theme
  scale_x_continuous(breaks = c(0, 1000, 2000)) +  # Customize breaks to show only 0, 1000, and 2000
  theme(axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

library(gridExtra)
# Combine the plots into one figure
combined_plot <- grid.arrange(p1, p2, ncol = 1, heights = c(1, 1))

# Save the combined plot as a TIFF file with the specified dimensions
ggsave("Figure2.tiff", plot = combined_plot, width = 170, height = 220, dpi = 300, units = "mm")