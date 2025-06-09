#G-compuatation function
#function for computing SPCE
#with U
Gcompute_spce_withU <- function(data) {
  # Fit seperated weibull AFT models for A=1 and A=0 by using sub-dataset
  model1 <- survreg(Surv(M, status) ~ x1 + x2 + U, data = subset(data, A == 1), dist = "weibull")
  model0 <- survreg(Surv(M, status) ~ x1 + x2 + U, data = subset(data, A == 0), dist = "weibull")
  
  # Predict survival probability at specific time point over entire dataset
  mu1 <- predict(model1, newdata = data, type = "lp")
  mu0 <- predict(model0, newdata = data, type = "lp")
  sigma1 <-  model1$scale
  sigma0 <-  model0$scale
  # S1 <- exp(-exp((log(t_pred)-mu1)/sigma1))
  # S0 <- exp(-exp((log(t_pred)-mu0)/sigma0))
  S1 <- exp(-(exp(-mu1) * t_pred)^(1/sigma1))
  S0 <- exp(-(exp(-mu0) * t_pred)^(1/sigma0))
  
  #log-normal aft
  # S1 <- 1 - pnorm((log(t_pred) - mu1) / model1$scale, mean=0, sd=1) #this only works for log-normal aft, not weibull
  # S0 <- 1 - pnorm((log(t_pred) - mu0) / model0$scale, mean=0, sd=1) 
  #S(t)= 1- phi((log(t)-mu)/sigma), phi is CDF of standard normal
  #OR S <- 1 - plnorm(t_pred, meanlog = mu, sdlog = model$scale) 
  
  return(mean(S1) - mean(S0)) # Estimate SPCE
}

Gcompute_spce_withoutU <- function(data) {
  # Fit seperated weibull AFT models for A=1 and A=0 by using sub-dataset
  model1 <- survreg(Surv(M, status) ~ x1 + x2 , data = subset(data, A == 1), dist = "weibull")
  model0 <- survreg(Surv(M, status) ~ x1 + x2 , data = subset(data, A == 0), dist = "weibull")
  # Predict survival probability at specific time point over entire dataset
  mu1 <- predict(model1, newdata = data, type = "lp")
  mu0 <- predict(model0, newdata = data, type = "lp")
  sigma1 <-  model1$scale
  sigma0 <-  model0$scale
  # S1 <- exp(-exp((log(t_pred)-mu1)/sigma1)) 
  # S0 <- exp(-exp((log(t_pred)-mu0)/sigma0))
  S1 <- exp(-(exp(-mu1) * t_pred)^(1/sigma1))
  S0 <- exp(-(exp(-mu0) * t_pred)^(1/sigma0))
  return(mean(S1) - mean(S0))
}