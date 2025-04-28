surv_stoEM_ipw <- function(data, zetat, zetaz, B = 100, theta = 0.5){
  t = data$M
  d = data$status
  Z = data$A
  X = data.matrix(data[,c("x1","x2")])
  nx = 2
  n = length(t)
  

  
  #Record coefficients with simulated U
  Coeft1 = Coeft1.se =spce = numeric(B)
  partialR2z = c()
  partialR2t1 = c()
  
  for (j in 1:B){
    Usim = SimulateU_surv(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset=FALSE)
    
    Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
    # Calculate partial R-sq of z ~ u | x
    Z.fit_reduced = glm(Z ~ X, family=binomial(link="probit"))
    # if (zetaz >= 0)
    #   partialR2z[j] = 1 - exp((Z.fit$deviance-Z.fit_reduced$deviance)/n)
    # else
    #   partialR2z[j] = - (1 - exp((Z.fit$deviance-Z.fit_reduced$deviance)/n))
    
    partialR2z[j] <- ifelse(zetaz >= 0,  1 - exp((Z.fit$deviance-Z.fit_reduced$deviance)/n), - (1 - exp((Z.fit$deviance-Z.fit_reduced$deviance)/n)))
    
    #ps = pnorm(cbind(1,X,Usim$U) %*% c(Z.fit$coefficients, zetaz))
    ps = Z.fit$fitted.values
    ipw = (sum(Z)/n) * Z/ps + (1-(sum(Z)/n))*(1-Z)/(1-ps)
    ipw = pmin(ipw, 10)
    ipw = pmax(ipw, 0.1)
    
    t1.ipw = coxph(Surv(t, d) ~ Z, weights = ipw, robust = TRUE)
    # Coeft1[j] = t1.ipw$coefficients
    # Coeft1.se[j] = summary(t1.ipw)$coefficients[,'robust se']
    spce[j] = mean(predict(t1.ipw, newdata=data.frame(Z=c(rep(1,n))),type = "survival", times = 2))-
                       mean(predict(t1.ipw, newdata=data.frame(Z=c(rep(0,n))),type = "survival", times = 2))
    
    
    # Calculate partial R-sq of (t, d) ~ u | x, z
    t1.fit = coxph(Surv(t, d) ~ X + Z + offset(zetat * Usim$U))
    t1.fit_reduced = coxph(Surv(t, d) ~ X + Z)
    logtest <- -2 * (t1.fit_reduced$loglik[2] - t1.fit$loglik[2])
    if (zetat >= 0)
      partialR2t1[j] = (1 - exp(-logtest/t1.fit$nevent))
    else
      partialR2t1[j] = - (1 - exp(-logtest/t1.fit$nevent))
  }
  
  # names(Coeft1) =  names(t1.ipw$coefficients)
  
  tau1 = mean(Coeft1)
  spce = mean(spce)
  # tau1.se = sqrt(mean((Coeft1.se)^2) + (1+1/B) * var(Coeft1))
  pR2z = mean(partialR2z)
  pR2t1 = mean(partialR2t1)
  
  return (list(tau1 = tau1, spce = spce, pR2z = pR2z, pR2t1 = pR2t1))
}