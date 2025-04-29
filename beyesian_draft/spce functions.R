#function for computing SPCE
compute_spce_withU <- function(data) {
  # Fit seperated weibull AFT models for A=1 and A=0 by using sub-dataset
  model1 <- survreg(Surv(M, status) ~ x1 + x2 + U, data = subset(data, A == 1), dist = "weibull")
  model0 <- survreg(Surv(M, status) ~ x1 + x2 + U, data = subset(data, A == 0), dist = "weibull")
  
  # Predict survival probability at specific time point over entire dataset
  #in survreg, T~weibull(lamda(scale) = exp(mu), shape = v = 1/sigma)
  #S(t)= exp(-(t/lamda)^v)
  mu1 <- predict(model1, newdata = data, type = "lp")
  mu0 <- predict(model0, newdata = data, type = "lp")
  sigma1 <-  model1$scale
  sigma0 <-  model0$scale
  scale1 <- exp(mu1)
  scale0 <- exp(mu0)
  shape1 <- 1 / sigma1
  shape0 <- 1 / sigma0
  
  S1 <- exp( - (t_pred / scale1)^shape1 )
  S0 <- exp( - (t_pred / scale0)^shape0 )
  
  
  return(mean(S1) - mean(S0)) # Estimate SPCE
}


compute_spce_withoutU <- function(data) {
  # Fit seperated weibull AFT models for A=1 and A=0 by using sub-dataset
  model1 <- survreg(Surv(M, status) ~ x1 + x2, data = subset(data, A == 1), dist = "weibull")
  model0 <- survreg(Surv(M, status) ~ x1 + x2, data = subset(data, A == 0), dist = "weibull")
  
  # Predict survival probability at specific time point over entire dataset
  mu1 <- predict(model1, newdata = data, type = "lp")
  mu0 <- predict(model0, newdata = data, type = "lp")
  sigma1 <-  model1$scale
  sigma0 <-  model0$scale
  scale1 <- exp(mu1)
  scale0 <- exp(mu0)
  shape1 <- 1 / sigma1
  shape0 <- 1 / sigma0
  
  S1 <- exp( - (t_pred / scale1)^shape1 )
  S0 <- exp( - (t_pred / scale0)^shape0 )
  
  return(mean(S1) - mean(S0)) # Estimate SPCE
}