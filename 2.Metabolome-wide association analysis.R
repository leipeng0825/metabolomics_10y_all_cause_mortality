#1-Load necessary RData file
load("~/Allcause_Mortality/allcausemortality.RData") 

# Read the UK Biobank CSV file into a dataframe
ukb_all_new <- read.csv("ukb_fulldata_final.csv")

# Check the columns in the dataset
column_names <- names(ukb_all_new)
print(column_names)

# 2-Random forest imputation for missing data
# Load the missRanger package for imputation
library(missRanger)

# Select the columns to be imputed
columns_to_impute <- c(1:12, 13:261)  #1:12 means 12 covariates, 13:261 means 249 metabolomic biomarkers

# Extract the subset of data that needs to be imputed
subset_to_impute <- ukb_all_new[, columns_to_impute]

# Perform random forest imputation on the selected columns using 100 trees
imputed_subset <- missRanger(subset_to_impute, num.trees = 100, seed = 1234)

# Create a new dataframe to store the imputed results
ukb_all_new_imputed <- ukb_all_new

# Place the imputed columns back into the original dataframe
ukb_all_new_imputed[, columns_to_impute] <- imputed_subset

# Check if there are any missing values in the imputed dataframe
summary(ukb_all_new_imputed)

# View the first few rows of the imputed dataframe
head(ukb_all_new_imputed)

# Check the columns in the imputed dataframe
column_names <- names(ukb_all_new_imputed)
print(column_names)

# Load the dplyr package for data manipulation
library(dplyr)

# 3-Data transformation
# Log-transform variables in columns 13 to 261
ukb_all_new_imputed_logall <- ukb_all_new_imputed %>%
  mutate(across(13:261, ~ log(.x)))

# Define z-score transformation function for a single column
z_score <- function(metab) {
  avg <- mean(metab, na.rm = TRUE)
  std <- sd(metab, na.rm = TRUE)
  
  # Perform z-score standardization
  metab <- (metab - avg) / std
  return(metab)
}

# Define z-score transformation function for a dataframe
z_tranb <- function(dataset, metab_col_start, metab_col_end) {
  # Create a copy of the dataset to store the standardized values
  z.data <- dataset
  
  # Loop through the specified columns and apply z-score transformation
  for (j in metab_col_start:metab_col_end) {
    z.data[, j] <- z_score(dataset[, j])
  }
  
  return(z.data)
}

# Perform z-score normalization on columns 13 to 261 in the dataframe ukb_all_new_imputed_logall
ukb_all_new_imputed_logall_z <- z_tranb(ukb_all_new_imputed_logall, 13, 261)

# View the first few rows of the dataframe after log transformation and z-score normalization
head(ukb_all_new_imputed_logall_z)

# Using the same analysis methods as applied to the UK Biobank database, we obtained a comparable dataset in the external validation with the ESTHER database: es_all_new_imputed_logall_z

# 4-Cox regression in total study population
# Load necessary packages
library(broom)
library(dplyr)
library(survival)

# Identify covariates, time-to-event variable, outcome variable, and metabolite columns
covariates <- colnames(ukb_all_new_imputed_logall_z)[c(1:12)]
time_var <- colnames(ukb_all_new_imputed_logall_z)[262]
outcome_var <- colnames(ukb_all_new_imputed_logall_z)[263]
metabolite_columns <- colnames(ukb_all_new_imputed_logall_z)[13:261]

# List to store Cox regression results
cox_results <- list()

# Run Cox regression model and add Schoenfeld test results
for (metab_col in metabolite_columns) {
  # Construct Cox model formula
  formula <- as.formula(paste("Surv(", time_var, ", ", outcome_var, ") ~ ", metab_col, " + ", paste(covariates, collapse = " + ")))
  
  # Fit Cox model
  cox_model <- coxph(formula, data = ukb_all_new_imputed_logall_z)
  
  # Get tidy summary of Cox model
  cox_summary <- tidy(cox_model, conf.int = TRUE)
  cox_summary <- cox_summary %>% filter(term == metab_col)  # Keep only metabolite-related result
  cox_summary$metabolite <- metab_col  # Add metabolite column
  
  # Perform Schoenfeld residual test
  schoenfeld_test <- cox.zph(cox_model)
  
  # Check if the metabolite is in the Schoenfeld test table
  if (metab_col %in% rownames(schoenfeld_test$table)) {
    schoenfeld_p_value <- schoenfeld_test$table[metab_col, "p"]
  } else {
    schoenfeld_p_value <- NA  # Handle cases where Schoenfeld test result is not available
  }
  
  # Add Schoenfeld test p-value to Cox summary
  cox_summary$schoenfeld_p_value <- schoenfeld_p_value
  
  # Store result
  cox_results[[metab_col]] <- cox_summary
}

# Combine the results into a data frame
results_df <- bind_rows(cox_results)

# Calculate RR and 95% confidence intervals
results_df <- results_df %>%
  mutate(RR = exp(estimate),
         lower_CI = exp(conf.low),
         upper_CI = exp(conf.high))

# Perform FDR adjustment
results_df$p.adjust <- p.adjust(results_df$p.value, method = "fdr")
results_df$log10_FDR_p.value <- -log10(results_df$p.adjust)
results_df$log2_RR <- log2(results_df$RR)

# Rename columns for clarity
results_df <- results_df %>%
  select(metabolite, estimate, std.error, RR, lower_CI, upper_CI, log2_RR, p.value, p.adjust, log10_FDR_p.value, schoenfeld_p_value)

# Save results to Excel
write.xlsx(results_df, "cox_ukb_all_249results.xlsx", rowNames = FALSE)

# 5-Create volcano plot for 249 metabolic biomarkers in the UK Biobank
# Load necessary packages
library(ggplot2)
library(ggrepel)
library(dplyr)

# Assume results_df is the dataframe containing Cox regression results
# Create logRR and log10_FDR.p.value columns
results_df <- results_df %>%
  mutate(logRR = log(RR),
         log10_FDR.p.value = -log10(p.adjust))

# Identify top 10 metabolites with RR greater than 1 and RR less than 1
top_metabolites_greater_than_1 <- results_df %>%
  filter(RR > 1) %>%
  arrange(desc(RR)) %>%
  head(10)

top_metabolites_less_than_1 <- results_df %>%
  filter(RR < 1) %>%
  arrange(RR) %>%
  head(10)

# Label these metabolites for visualization
results_df$label <- ifelse(results_df$metabolite %in% c(top_metabolites_greater_than_1$metabolite, 
                                                        top_metabolites_less_than_1$metabolite), 
                           results_df$metabolite, "")

# Create a new column to categorize data points for coloring based on significance and effect size
results_df$Risk_Classification <- ifelse(results_df$p.adjust < 0.05, 
                                         ifelse(results_df$logRR > 0, "Positive", "Negative"), 
                                         "Nonsignificant")

# Count the number of points in each Risk_Classification category
label_summary <- table(results_df$Risk_Classification)
print(label_summary)

# Dynamically create labels for the color categories based on the summary
color_labels <- paste(names(label_summary), "(n=", as.numeric(label_summary), ")", sep = "")

# Create the volcano plot
volcano_plot1 <- ggplot(results_df, aes(x = logRR, y = log10_FDR.p.value, color = Risk_Classification)) +
  geom_point(shape = 19, alpha = 1) +  # Plot points
  geom_text_repel(aes(label = label), max.overlaps = 20, size = 3, alpha = 1) +  # Add labels for top metabolites
  geom_vline(xintercept = 0, linetype = 4, color = "black") +  # Vertical line at log(HR) = 0
  geom_hline(yintercept = -log10(0.05), linetype = 4, color = "red", size = 0.8) +  # Horizontal line at p = 0.05
  labs(title = "Associations between 249 metabolic biomarkers \n and 10-year all-cause mortality in UK Biobank (12,347 deceased)",
       x = "log(HR)",
       y = "-log10(FDR-adjusted p-value)") +
  scale_color_manual(values = c("#104E8B", "#A6A6A6", "#8B1A1A"), 
                     labels = color_labels) +  # Set color and labels for each category
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),  # Black border
        panel.background = element_rect(fill = "white"),  # White background
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.text = element_text(size = 10),  # Set axis text size
        panel.grid = element_blank(),  # Remove grid lines
        legend.position = c(0.13, 0.95),   # Position the legend in the upper left
        axis.ticks = element_line())  # Show axis ticks

# Display the volcano plot
print(volcano_plot1)

# Save the volcano plot to a PDF file
ggsave("cox_ukb_all_volcano_249biomarkers.pdf", plot = volcano_plot1, width = 8, height = 7)

# Extract significant 224 metabolic biomarkers based on FDR-adjusted p-value
significant_metabolites <- results_df %>%
  filter(p.adjust < 0.05) %>%
  pull(metabolite)

# Validate these 224 metabolic biomarkers in the ESTHER Study
# List to store the new analysis results
new_cox_results <- list()

# Run new Cox regression models and perform Schoenfeld residual test in the ESTHER Study
for (metab_col in significant_metabolites) {
  # Create the Cox model formula
  formula <- as.formula(paste("Surv(", time_var, ", ", outcome_var, ") ~ ", metab_col, " + ", paste(covariates, collapse = " + ")))
  
  # Fit the Cox model using ESTHER Study data
  new_cox_model <- coxph(formula, data = es_all_new_imputed_logall_z)  # ESTHER Study data
  
  # Summarize the model output and filter for metabolite result
  new_cox_summary <- tidy(new_cox_model, conf.int = TRUE)
  new_cox_summary <- new_cox_summary %>% filter(term == metab_col)  # Keep only metabolite-related results
  new_cox_summary$metabolite <- metab_col  # Add metabolite column
  
  # Perform Schoenfeld residual test
  schoenfeld_test <- cox.zph(new_cox_model)
  
  # Check if the metabolite is present in the Schoenfeld test result table
  if (metab_col %in% rownames(schoenfeld_test$table)) {
    schoenfeld_p_value <- schoenfeld_test$table[metab_col, "p"]
  } else {
    schoenfeld_p_value <- NA  # Handle cases where the Schoenfeld test result is not available
  }
  
  # Add Schoenfeld test p-value to the summary
  new_cox_summary$schoenfeld_p_value <- schoenfeld_p_value
  
  # Store the result in the list
  new_cox_results[[metab_col]] <- new_cox_summary
}

# Combine the new results into a dataframe
new_results_df <- bind_rows(new_cox_results)

# Calculate RR and 95% confidence intervals
new_results_df <- new_results_df %>%
  mutate(RR = exp(estimate),
         lower_CI = exp(conf.low),
         upper_CI = exp(conf.high))

# Perform FDR adjustment on p-values
new_results_df$p.adjust <- p.adjust(new_results_df$p.value, method = "fdr")
new_results_df$log10_FDR_p.value <- -log10(new_results_df$p.adjust)
new_results_df$log2_RR <- log(new_results_df$RR)

# Rename columns for better understanding
new_results_df <- new_results_df %>%
  select(metabolite, estimate, std.error, RR, lower_CI, upper_CI, log2_RR, p.value, p.adjust, log10_FDR_p.value, schoenfeld_p_value)

# Save the new results to an Excel file
write.xlsx(new_results_df, "cox_es_all_224results.xlsx", rowNames = FALSE)

# 6-Create volcano plot for 224 metabolic biomarkers in the ESTHER Study
# Load necessary packages
library(ggplot2)
library(ggrepel)
library(dplyr)

# Create logRR and log10_FDR.p.value columns
results_df <- new_results_df %>%
  mutate(logRR = log(RR),
         log10_FDR.p.value = -log10(p.adjust))

# Identify top metabolites with RR greater than 1 and RR less than 1
top_metabolites_greater_than_1 <- results_df %>%
  filter(RR > 1, p.adjust < 0.05) %>%
  arrange(desc(RR)) %>%
  head(10)

top_metabolites_less_than_1 <- results_df %>%
  filter(RR < 1, p.adjust < 0.05) %>%
  arrange(RR) %>%
  head(10)

# Label these metabolites
results_df$label <- ifelse(results_df$metabolite %in% c(top_metabolites_greater_than_1$metabolite, 
                                                        top_metabolites_less_than_1$metabolite), 
                           results_df$metabolite, "")

# Create a new column for coloring based on significance and effect size
results_df$Risk_Classification <- ifelse(results_df$log10_FDR.p.value > -log10(0.05), 
                                         ifelse(results_df$logRR > 0, "Positive", "Negative"), 
                                         "Nonsignificant")

# Count the number of points in each Risk_Classification category
label_summary <- table(results_df$Risk_Classification)
print(label_summary)

# Dynamically create color labels based on classification count
color_labels <- paste(names(label_summary), "(n=", as.numeric(label_summary), ")", sep = "")

# Plot the volcano plot
volcano_plot2 <- ggplot(results_df, aes(x = logRR, y = log10_FDR.p.value, color = Risk_Classification)) +
  geom_point(shape = 19, alpha = 1) +
  geom_text_repel(aes(label = label), max.overlaps = 20, size = 3, alpha = 1) +
  geom_vline(xintercept = 0, linetype = 4, color = "black") +  # Vertical reference line at log(HR) = 0
  geom_hline(yintercept = -log10(0.05), linetype = 4, color = "red", size = 0.8) +  # Horizontal reference line at p = 0.05
  labs(title = "Associations between 224 metabolic biomarkers \n and 10-year all-cause mortality in ESTHER Study (804 deceased)",
       x = "log(HR)",
       y = "-log10(FDR-adjusted p-value)") +
  scale_color_manual(values = c("#104E8B", "#A6A6A6", "#8B1A1A"), 
                     labels = color_labels) +  # Dynamically set color labels
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),  # Black border
        panel.background = element_rect(fill = "transparent"),  # Transparent background
        plot.title = element_text(hjust = 0.5),  # Center title
        axis.text = element_text(size = 10),  # Set axis text size
        panel.grid = element_blank(), # Remove grid lines
        legend.position = c(0.13, 0.95),   # Position legend in the upper left corner
        axis.ticks = element_line())  # Show axis ticks

# Display the volcano plot
print(volcano_plot2)

# Save the plot to a PDF file
ggsave("cox_es_all_volcano_224biomarkers.pdf", plot = volcano_plot2, width = 8, height = 7)
