# 1-Load necessary RData file
load("~/Allcause_Mortality/allcausemortality.RData") 

# Set seed for reproducibility
set.seed(234)

# Generate random indices for splitting data into training and validation sets
train_index <- sample(seq_len(nrow(ukb_all_new_imputed_logall_z)), size = 0.7 * nrow(ukb_all_new_imputed_logall_z))

# 2-Split data into training set and validation set
ukb_train_set <- ukb_all_new_imputed_logall_z[train_index, ]
ukb_validation_set <- ukb_all_new_imputed_logall_z[-train_index, ]

# Print first few rows of training and validation sets
print(head(ukb_train_set))
print(head(ukb_validation_set))

# 3-Extract subgroups
# Extract 4 subgroups from training set
ukb_train_male5059 <- subset(ukb_train_set, sex == 1 & (age >= 50 & age <= 59))
ukb_train_male6069 <- subset(ukb_train_set, sex == 1 & (age >= 60 & age <= 69))
ukb_train_female5059 <- subset(ukb_train_set, sex == 0 & (age >= 50 & age <= 59))
ukb_train_female6069 <- subset(ukb_train_set, sex == 0 & (age >= 60 & age <= 69))

# Extract 4 subgroups from internal validation set
ukb_validation_male5059 <- subset(ukb_validation_set, sex == 1 & (age >= 50 & age <= 59))
ukb_validation_male6069 <- subset(ukb_validation_set, sex == 1 & (age >= 60 & age <= 69))
ukb_validation_female5059 <- subset(ukb_validation_set, sex == 0 & (age >= 50 & age <= 59))
ukb_validation_female6069 <- subset(ukb_validation_set, sex == 0 & (age >= 60 & age <= 69))

# Extract 4 subgroups from external validation set
es_all_new_male5059 <- subset(es_all_new_imputed_logall_z, sex == 2 & (age >= 50 & age <= 59))
es_all_new_male6069 <- subset(es_all_new_imputed_logall_z, sex == 2 & (age >= 60 & age <= 69))
es_all_new_female5059 <- subset(es_all_new_imputed_logall_z, sex == 1 & (age >= 50 & age <= 59))
es_all_new_female6069 <- subset(es_all_new_imputed_logall_z, sex == 1 & (age >= 60 & age <= 69))

# Check the number of rows for each subgroup
print(nrow(ukb_train_male5059))
print(nrow(ukb_train_male6069))
print(nrow(ukb_train_female5059))
print(nrow(ukb_train_female6069))

print(nrow(ukb_validation_male5059))
print(nrow(ukb_validation_male6069))
print(nrow(ukb_validation_female5059))
print(nrow(ukb_validation_female6069))

print(nrow(es_all_new_male5059))
print(nrow(es_all_new_male6069))
print(nrow(es_all_new_female5059))
print(nrow(es_all_new_female6069))

# 4-Define a function to perform Bootstrap LASSO analysis for a given dataset
perform_bootstrap_lasso <- function(data, group_name) {
  # Set the dependent variable (y)
  data$followup_year <- round(data$followup_year)
  y <- Surv(data$followup_year, data$mortality10y)
  
  # Set the independent variables (x)
  x <- data[, 13:261]  # Columns 13 to 261 are the metabolic biomarkers
  
  # Combine independent variables and dependent variables into a data frame
  analysis_data <- as.data.frame(cbind(x, followup_year = data$followup_year, mortality10y = data$mortality10y))
  
  # Define the Bootstrap LASSO function
  bootstrap_lasso <- function(data, indices) {
    d <- data[indices, ]
    x_boot <- as.matrix(d[, -c(ncol(d)-1, ncol(d))])  # Exclude followup_year and mortality10y
    y_boot <- Surv(d$followup_year, d$mortality10y)
    
    cv_lasso <- cv.glmnet(x_boot, y_boot, alpha = 1, family = "cox", nfolds = 10)
    best_lambda <- cv_lasso$lambda.1se
    lasso_model <- glmnet(x_boot, y_boot, alpha = 1, family = "cox", lambda = best_lambda)
    
    return(as.vector(coef(lasso_model)) != 0)
  }
  
  # Perform 1000 Bootstrap iterations
  library(boot)
  set.seed(234)
  results <- boot(analysis_data, statistic = bootstrap_lasso, R = 1000)
  
  # Count how many times each metabolite was selected
  selected_counts <- apply(results$t, 2, sum)
  
  # Select metabolites that were chosen in at least 95% of the Bootstrap samples
  selected_metabolites <- colnames(x)[selected_counts >= 0.95 * 1000]
  
  # Print the selected metabolites
  print(paste("Selected metabolites for", group_name, "(chosen in at least 95% of the bootstrap samples):"))
  print(selected_metabolites)
  
  # Save the selected metabolites to a CSV file
  output_filename <- paste0("bootstrap_selected_metabolites_", group_name, ".csv")
  write.csv(selected_metabolites, output_filename, row.names = FALSE)
}

# Apply the function to all four subgroups
# For male aged 50-59
perform_bootstrap_lasso(ukb_train_male5059, "male5059")
# For male aged 60-69
perform_bootstrap_lasso(ukb_train_male6069, "male6069")
# For female aged 50-59
perform_bootstrap_lasso(ukb_train_female5059, "female5059")
# For female aged 60-69
perform_bootstrap_lasso(ukb_train_female6069, "female6069")


# 5-Cox regression analysis and forest plots
library(survival)
library(ggplot2)
library(dplyr)

# Function to perform multivariable Cox regression
cox_multivariable_model <- function(data, biomarkers, covariates, outcome, followup) {
  data <- data %>% mutate(Surv_obj = Surv(get(followup), get(outcome)))
  formula <- as.formula(paste("Surv_obj ~", paste(c(biomarkers, covariates), collapse = " + ")))
  cox_model <- coxph(formula, data = data)
  summary_cox <- summary(cox_model)

  result <- data.frame(
    Variable = rownames(summary_cox$coefficients),
    HR = summary_cox$coefficients[, "exp(coef)"],
    CI_Lower = summary_cox$conf.int[, "lower .95"],
    CI_Upper = summary_cox$conf.int[, "upper .95"],
    P_Value = summary_cox$coefficients[, "Pr(>|z|)"]
  )
  return(result)
}

# Function to combine results across datasets and export to CSV
combine_and_export_full <- function(results_train, results_validation, results_external, filename) {
  results_train$Dataset <- "Training"
  results_validation$Dataset <- "Validation"
  results_external$Dataset <- "External"
  all_results <- rbind(results_train, results_validation, results_external)
  write.csv(all_results, file = filename, row.names = FALSE)
}

# Function to plot forest plots for biomarkers from multivariable Cox model
plot_forest_from_full_model <- function(filename, target_vars, x_labels, outcome_name, tiff_filename = NULL) {
  results <- read.csv(filename)
  results <- results[results$Variable %in% target_vars, ]  # Keep only biomarkers for plotting

  results <- results %>%
    mutate(Variable = factor(Variable, levels = target_vars),
           Dataset = factor(Dataset, levels = c("Training", "Validation", "External"))) %>%
    arrange(Variable, Dataset)

  p <- ggplot(results, aes(x = Variable, y = HR, ymin = CI_Lower, ymax = CI_Upper, color = Dataset)) +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.15) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_continuous(limits = c(0.2, 2.5), breaks = seq(0.2, 2.5, by = 0.2)) +
    scale_x_discrete(labels = x_labels) +
    labs(title = outcome_name, x = "", y = "HR per 1-SD increment") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    scale_color_manual(values = c("Training" = "#000000", "Validation" = "#CD3333", "External" = "#1874CD"),
                       labels = c("UK Biobank (70%)", "UK Biobank (30%)", "ESTHER Study"))

  if (!is.null(tiff_filename)) {
    ggsave(filename = tiff_filename, plot = p, device = "tiff", width = 7, height = 5, units = "in", dpi = 600)
  } else {
    print(p)
  }
}

# 6-Define covariates, outcome variables, and follow-up time
covariates12 <- c("age", "sex", "education", "smoking", "activity", "alcohol",
                  "BMI", "hypertension", "diabetes", "dyslipidemia", "CVD", "cancer")
outcome <- c("mortality10y")
followup_year <- "followup_year"

# Function to check required columns exist
check_columns <- function(data, outcomes, followup) {
  all_columns <- c(outcomes, followup)
  if (!all(all_columns %in% colnames(data))) {
    stop("Required outcome or follow-up columns are missing from the dataset")
  }
}

# Loop through each subgroup and run multivariable Cox models
for (group in names(group_info)) {
  info <- group_info[[group]]

  # Validate column presence
  check_columns(info$train_data, outcome, followup_year)
  check_columns(info$validation_data, outcome, followup_year)
  check_columns(info$external_data, outcome, followup_year)

  # Perform multivariable Cox analysis
  for (i in 1:length(outcome)) {
    outcome_var <- outcome[i]
    result_train <- cox_multivariable_model(info$train_data, info$biomarkers, covariates12, outcome_var, followup_year)
    result_valid <- cox_multivariable_model(info$validation_data, info$biomarkers, covariates12, outcome_var, followup_year)
    result_external <- cox_multivariable_model(info$external_data, info$biomarkers, covariates12, outcome_var, followup_year)

    filename <- paste0("cox_full_model_results_", group, "_", i, ".csv")
    combine_and_export_full(result_train, result_valid, result_external, filename)

    # Plot forest plot (biomarkers only)
    plot_forest_from_full_model(
      filename = filename,
      target_vars = info$biomarkers,
      x_labels = info$x_labels,
      outcome_name = paste0("All-cause mortality (", group, ")"),
      tiff_filename = paste0("forestplot_fullmodel_", group, ".tiff")
    )
  }
}
