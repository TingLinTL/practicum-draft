compute_ipw_spce_withU <- function(data, time_var, event_var, treat_var, t_eval) {
  
  #treatment variable
  A <- data[[treat_var]]
  covariates <- c("x1", "x2", "U")
    
  # Estimate propensity scores
  ps_formula <- as.formula(paste(treat_var, "~", paste(covariates, collapse = " + ")))
  ps_model <- glm(ps_formula, data = data, family = binomial)
  ps <- predict(ps_model, type = "response")
  
  #stabilized weights
  pA <- mean(A)
  sw <- ifelse(A == 1, pA / ps, (1 - pA) / (1 - ps))
  data$sw <- sw
  
  #weighted Kaplan-Meier estimator
  surv_obj <- Surv(time = data[[time_var]], event = data[[event_var]])
  formula_km <- as.formula(paste("surv_obj ~", treat_var))
  km_fit <- survfit(formula_km, data = data, weights = sw)
  
  #survival estimates at specified time
  km_summary <- summary(km_fit, times = t_eval)
  surv_probs <- km_summary$surv
  group_labels <- sub(".*=", "", km_summary$strata)
  
  # Extract survival estimates for A = 1 and A = 0
  S1 <- as.numeric(surv_probs[group_labels == "1"])
  S0 <- as.numeric(surv_probs[group_labels == "0"])
  
  # Compute SPCE
  SPCE <- S1 - S0
  
  # Return results
  return(list(
    S1 = S1,
    S0 = S0,
    SPCE = SPCE
  ))
}


compute_ipw_spce_withoutU <- function(data, time_var, event_var, treat_var, t_eval) {
  
  #treatment variable
  A <- data[[treat_var]]
  covariates <- c("x1", "x2")
  
  # Estimate propensity scores
  ps_formula <- as.formula(paste(treat_var, "~", paste(covariates, collapse = " + ")))
  ps_model <- glm(ps_formula, data = data, family = binomial)
  ps <- predict(ps_model, type = "response")
  
  #stabilized weights
  pA <- mean(A)
  sw <- ifelse(A == 1, pA / ps, (1 - pA) / (1 - ps))
  data$sw <- sw
  
  #weighted Kaplan-Meier estimator
  surv_obj <- Surv(time = data[[time_var]], event = data[[event_var]])
  formula_km <- as.formula(paste("surv_obj ~", treat_var))
  km_fit <- survfit(formula_km, data = data, weights = sw)
  
  #survival estimates at specified time
  km_summary <- summary(km_fit, times = t_eval)
  surv_probs <- km_summary$surv
  group_labels <- sub(".*=", "", km_summary$strata)
  
  # Extract survival estimates for A = 1 and A = 0
  S1 <- as.numeric(surv_probs[group_labels == "1"])
  S0 <- as.numeric(surv_probs[group_labels == "0"])
  
  # Compute SPCE
  SPCE <- S1 - S0
  
  # Return results
  return(list(
    S1 = S1,
    S0 = S0,
    SPCE = SPCE
  ))
}