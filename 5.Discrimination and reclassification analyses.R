# The datasets of main subgroups in this study include:
#  ukb_validation_male5059
#  ukb_validation_male6069
#  ukb_validation_female5059
#  ukb_validation_female6069
#  es_all_new_male5059
#  es_all_new_male6069
#  es_all_new_female5059
#  es_all_new_female6069

# Load necessary RData file
load("~/Allcause_Mortality/allcausemortality.RData") 

# 1- Load Required Packages
packages <- c("survival", "pROC", "ggplot2", "nricens")
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

# 2- C-statistic and CI Calculation
calculate_model_aucs_with_ci <- function(data, covariate_indices, metabolite_indices, outcome_col, time_col) {
  covariates <- names(data)[covariate_indices]
  metabolites <- names(data)[metabolite_indices]

  formula_cov <- as.formula(paste("Surv(", time_col, ",", outcome_col, ") ~", paste(covariates, collapse = " + ")))
  formula_metab <- as.formula(paste("Surv(", time_col, ",", outcome_col, ") ~", paste(metabolites, collapse = " + ")))
  formula_full <- as.formula(paste("Surv(", time_col, ",", outcome_col, ") ~", paste(c(covariates, metabolites), collapse = " + ")))

  model_cov <- coxph(formula_cov, data = data)
  model_metab <- coxph(formula_metab, data = data)
  model_full <- coxph(formula_full, data = data)

  lp_cov <- predict(model_cov, newdata = data, type = "lp")
  lp_metab <- predict(model_metab, newdata = data, type = "lp")
  lp_full <- predict(model_full, newdata = data, type = "lp")

  roc_cov <- roc(data[[outcome_col]], lp_cov)
  roc_metab <- roc(data[[outcome_col]], lp_metab)
  roc_full <- roc(data[[outcome_col]], lp_full)

  cstat_cov <- auc(roc_cov)
  cstat_metab <- auc(roc_metab)
  cstat_full <- auc(roc_full)

  ci_cov <- ci.auc(roc_cov)
  ci_metab <- ci.auc(roc_metab)
  ci_full <- ci.auc(roc_full)

  roc_data <- data.frame(
    TPR = c(auc_cov, auc_metab, auc_full),
    FPR = c(1 - auc_cov, 1 - auc_metab, 1 - auc_full),
    Model = factor(c("Covariates-only Model, C-statistic", "Metabolites-only Model, C-statistic", "Full Model, C-statistic"))
  )

  ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(size = 0.5) +
    geom_abline(linetype = "dashed", color = "black") +
    labs(title = "", x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal() +
    scale_color_manual(values = c("blue", "green", "red")) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      legend.position = c(0.57, 0.10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    )

  return(list(
    cstat_covariates = cstat_cov,
    ci_covariates = ci_cov,
    cstat_metabolites = cstat_metab,
    ci_metabolites = ci_metab,
    cstat_full = cstat_full,
    ci_full = ci_full
  ))
}

# 3- IDI and NRI Calculation
calculate_idi_nri <- function(data, covariate_indices, metabolite_indices, outcome_col, time_col, cutoffs = c(0, 0.05, 0.10, 1)) {
  covariates <- data[, covariate_indices]
  metabolites <- data[, metabolite_indices]
  all_predictors <- data[, c(covariate_indices, metabolite_indices)]

  data[[time_col]] <- as.numeric(as.character(data[[time_col]]))
  data[[outcome_col]] <- as.numeric(as.character(data[[outcome_col]]))

  model_base <- coxph(Surv(data[[time_col]], data[[outcome_col]]) ~ ., data = covariates)
  model_full <- coxph(Surv(data[[time_col]], data[[outcome_col]]) ~ ., data = all_predictors)

  lp_base <- predict(model_base, newdata = data, type = "lp")
  lp_full <- predict(model_full, newdata = data, type = "lp")

  risk_base <- 1 - exp(-lp_base)
  risk_full <- 1 - exp(-lp_full)

  rec_data <- data.frame(event = data[[outcome_col]], risk_base, risk_full)

  result <- reclassification(
    data = rec_data,
    cOutcome = 1,
    predrisk1 = risk_base,
    predrisk2 = risk_full,
    cutoff = cutoffs,
    updown = "category"
  )

  result_continuous <- reclassification(
    data = rec_data,
    cOutcome = 1,
    predrisk1 = risk_base,
    predrisk2 = risk_full,
    cutoff = 0,
    updown = "diff"
  )

  return(list(
    categorical_nri_idi = result,
    continuous_nri_idi = result_continuous
  ))
}

# 4- Validation: for example, ukb_validation_male5059.
covariate_indices <- c(1:13)
metabolite_indices <- c(14:33)
outcome_col <- "mortality10y"
time_col <- "followup_year"

data_example <- ukb_validation_male5059

auc_results <- calculate_model_aucs_with_ci(
  data_example,
  covariate_indices,
  metabolite_indices,
  outcome_col,
  time_col
)

auc_summary <- data.frame(
  Model = c("Covariates", "Metabolites", "Full"),
  C_statistic = c(auc_results$cstat_covariates, auc_results$cstat_metabolites, auc_results$cstat_full),
  CI_Lower = c(auc_results$ci_covariates[1], auc_results$ci_metabolites[1], auc_results$ci_full[1]),
  CI_Upper = c(auc_results$ci_covariates[3], auc_results$ci_metabolites[3], auc_results$ci_full[3])
)
write.csv(auc_summary, "auc_comparison_results.csv", row.names = FALSE)

nri_idi_results <- calculate_idi_nri(
  data_example,
  covariate_indices,
  metabolite_indices,
  outcome_col,
  time_col
)

capture.output(nri_idi_results$categorical_nri_idi, file = "categorical_nri_idi.txt")
capture.output(nri_idi_results$continuous_nri_idi, file = "continuous_nri_idi.txt")
