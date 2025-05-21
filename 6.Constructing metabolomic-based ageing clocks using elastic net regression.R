# Load necessary RData file
load("~/Allcause_Mortality/allcausemortality.RData") 

# 1-Complete R code for constructing MetaboMR clock1 & clock2 and Cox evaluation across subgroups

library(glmnet)
library(survival)
library(dplyr)

# 2-Define covariates and outcome
covariates12 <- c("age", "sex", "education", "smoking", "activity", "alcohol",
                  "BMI", "hypertension", "diabetes", "dyslipidemia", "CVD", "cancer")
outcome <- "mortality10y"
followup_year <- "followup_year"

# 3-List of subgroups and their selected biomarkers
subgroup_info <- list(
  male5059 = list(train = ukb_train_male5059, valid = ukb_validation_male5059, external = es_all_new_male5059,
                  biomarkers = c("Acetate","Acetone","Albumin","bOHbutyrate","Citrate","Creatinine","Glucose","GlycA",
                                "HDL_size","His","IDL_CE_pct","L_LDL_CE_pct", "LA_pct","Lactate", "Omega_6_by_Omega_3",
                                "Tyr", "Val","VLDL_size","XL_HDL_FC", "XXL_VLDL_PL_pct")),
  male6069 = list(train = ukb_train_male6069, valid = ukb_validation_male6069, external = es_all_new_male6069,
                  biomarkers = c("Acetate","Acetone","Albumin","bOHbutyrate", "Citrate","Creatinine","Glucose","GlycA",
                                "His","IDL_CE_pct","LA_pct","Omega_6_by_Omega_3","S_HDL_CE","S_LDL_CE", "Tyr","Val",
                                "VLDL_size","XL_HDL_FC")),
  female5059 = list(train = ukb_train_female5059, valid = ukb_validation_female5059, external = es_all_new_female5059,
                    biomarkers = c("Acetoacetate", "Albumin", "Creatinine", "Glucose", "GlycA", "His", "LA_pct", "Leu",
                                  "M_LDL_TG_pct", "Omega_6_by_Omega_3", "Val", "VLDL_size")),
  female6069 = list(train = ukb_train_female6069, valid = ukb_validation_female6069, external = es_all_new_female6069,
                    biomarkers = c("Acetoacetate", "Albumin", "Creatinine", "Glucose", "GlycA", "His", "LA_pct",
                                  "M_LDL_TG_pct", "Omega_6_by_Omega_3", "PUFA", "S_HDL_CE", "Val", "VLDL_size"))
)

# Store MetAA values and results
all_metaa_results <- list()
cox_results_table <- data.frame()
cor_results <- data.frame()
metaa_summary <- data.frame()

for (group in names(subgroup_info)) {
  info <- subgroup_info[[group]]
  train <- info$train
  valid <- info$valid
  external <- info$external
  biomarkers <- info$biomarkers
  predictors_clock1 <- biomarkers
  predictors_clock2 <- c(biomarkers, covariates12)

  build_clock <- function(data, predictors, model_name) {
    data_sub <- data[, c(predictors, "age", outcome, followup_year)]
    data_sub <- na.omit(data_sub)
    X <- scale(data_sub[, predictors])
    y <- data_sub$age

    fit <- cv.glmnet(X, y, alpha = 0.5)
    raw_age <- predict(fit, newx = X, s = "lambda.min")

    calib <- lm(raw_age ~ y)
    calibrated_age <- predict(calib)
    metaa <- calibrated_age - y

    # Pearson correlation between clock and age
    cor_test <- cor.test(as.numeric(calibrated_age), y)
    cor_results <<- rbind(cor_results, data.frame(
      Group = group,
      Clock = model_name,
      r = round(cor_test$estimate, 3),
      p = signif(cor_test$p.value, 3)
    ))

    return(list(clock = as.numeric(calibrated_age), metaa = as.numeric(metaa), model = fit,
                predictors = predictors, calib = calib, center = attr(X, "scaled:center"),
                scale = attr(X, "scaled:scale")))
  }

  extract_cox_result <- function(model, metaa_label, dataset_label) {
    s <- summary(model)
    data.frame(
      Group = group,
      Dataset = dataset_label,
      Clock = metaa_label,
      HR = round(s$conf.int[1, "exp(coef)"], 3),
      CI_Lower = round(s$conf.int[1, "lower .95"], 3),
      CI_Upper = round(s$conf.int[1, "upper .95"], 3),
      P_value = signif(s$coefficients[1, "Pr(>|z|)"], 3)
    )
  }

  get_metaa_stats <- function(values, label, dataset) {
    data.frame(
      Group = group,
      Clock = label,
      Dataset = dataset,
      Mean = round(mean(values, na.rm = TRUE), 3),
      SD = round(sd(values, na.rm = TRUE), 3)
    )
  }

  clock1 <- build_clock(train, predictors_clock1, "MetaboMR_clock1")
  train$MetAA1 <- clock1$metaa
  train$MetaboMR_clock1 <- clock1$clock

  clock2 <- build_clock(train, predictors_clock2, "MetaboMR_clock2")
  train$MetAA2 <- clock2$metaa
  train$MetaboMR_clock2 <- clock2$clock

  all_metaa_results[[group]] <- list(train = train, clock1 = clock1, clock2 = clock2)

  m1 <- coxph(Surv(train[[followup_year]], train[[outcome]]) ~ MetAA1, data = train)
  m2 <- coxph(as.formula(paste("Surv(", followup_year, ",", outcome, ") ~ MetAA1 +", paste(covariates12, collapse = "+"))), data = train)
  m3 <- coxph(Surv(train[[followup_year]], train[[outcome]]) ~ MetAA2, data = train)

  cox_results_table <- rbind(cox_results_table,
                             extract_cox_result(m1, "MetAA1_crude", "Train"),
                             extract_cox_result(m2, "MetAA1_adjusted", "Train"),
                             extract_cox_result(m3, "MetAA2_crude", "Train"))

  metaa_summary <- rbind(metaa_summary,
                         get_metaa_stats(clock1$metaa, "MetAA1", "Train"),
                         get_metaa_stats(clock2$metaa, "MetAA2", "Train"))

  predict_and_eval <- function(test_data, clock, name, clock_label) {
    test_data_sub <- test_data[, c(clock$predictors, "age", outcome, followup_year)]
    test_data_sub <- na.omit(test_data_sub)
    X_test <- scale(test_data_sub[, clock$predictors], center = clock$center, scale = clock$scale)
    raw_pred <- predict(clock$model, newx = X_test, s = "lambda.min")
    calibrated <- predict(clock$calib, newdata = data.frame(y = test_data_sub$age))
    metaa <- as.numeric(calibrated - test_data_sub$age)

    cor_test <- cor.test(as.numeric(calibrated), test_data_sub$age)
    cor_results <<- rbind(cor_results, data.frame(
      Group = group,
      Clock = clock_label,
      r = round(cor_test$estimate, 3),
      p = signif(cor_test$p.value, 3),
      Dataset = name
    ))

    metaa_summary <<- rbind(metaa_summary, get_metaa_stats(metaa, clock_label, name))

    m <- coxph(Surv(test_data_sub[[followup_year]], test_data_sub[[outcome]]) ~ metaa, data = test_data_sub)
    cox_results_table <<- rbind(cox_results_table, extract_cox_result(m, clock_label, name))
  }

  predict_and_eval(valid, clock1, "Internal", "MetAA1")
  predict_and_eval(external, clock1, "External", "MetAA1")
  predict_and_eval(valid, clock2, "Internal", "MetAA2")
  predict_and_eval(external, clock2, "External", "MetAA2")
}

# 4-Export all results
write.csv(cox_results_table, "cox_metaa_all_results.csv", row.names = FALSE)
write.csv(cor_results, "clock_age_correlation_results.csv", row.names = FALSE)
write.csv(metaa_summary, "metaa_distribution_summary.csv", row.names = FALSE)
