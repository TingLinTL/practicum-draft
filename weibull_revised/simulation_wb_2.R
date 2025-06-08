#Simulation
library(survival)
set.seed(2000)
n_sim <-200 #number of simulation
n <- 1000 #sample size 
tau <- 5.5 #maximum follow-up time
t_pred <- 2  #time point for estimates 
true_spce <- 0.1901985 #true SPCE at time=2 from numeric iteration for weibull aft model
nboot<-1000
spce_withU <- spce_withoutU <- spce_huang <- numeric(n_sim)
spce_withU.b <- spce_withoutU.b <- spce_huang.b <- ipw_spce_withU.b <- ipw_spce_withoutU.b <- rep(NA, nboot)
spce_withU.sd <- spce_withoutU.sd <- spce_huang.sd <- spce_withU_iptw.sd <- spce_withoutU_iptw.sd <- numeric(n_sim)
cp_withU <- cp_withoutU <- cp_huang <- logical(n_sim)
ipw_spce_withU <- ipw_spce_withoutU <- numeric(n_sim)


#function for computing SPCE
#with U
compute_spce_withU <- function(data) {
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

compute_spce_withoutU <- function(data) {
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

#simulation
for (i in 1:n_sim) {
  #Generate X1, X2, U
  X1 <- rbinom(n, size = 1, prob = 0.4) #X1~Bern(0.4)
  X2 <- rbinom(n, size = 1, prob = 0.6) #X2~Bern(0.6)
  U <- rbinom(n, size = 1, prob = 0.5) #U~Bern(0.5), independent with X1 and X2
  
  # Treatment assignment, logit(P(A=1|X1,X2,U))
  alpha_0 <- 0.1
  alpha_x1 <- 0.3
  alpha_x2 <- 0.5
  alpha_u <- -0.8
  pi_A <- 1 / (1 + exp(-(alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U)))
  A <- rbinom(n, size = 1, prob = pi_A)
  
  #Potential outcomes
  #event
  eta_intercept <- 0.2 #A=0
  eta_intercept_a1<- 0.7 #A=1
  eta_x1 <- -0.1
  eta_x2 <- 0.4
  eta_u <- -0.8
  
  #C_a0, C_a1 potential censored time
  C_a0 <- C_a1 <- runif(n, 0.1, 5.5) #non-informative censoring
  
  #D_a0, D_a1 are genrated from weibull aft model
  #D_a0,D_a1 potential event time
  #inverse transformation method
  #t_i* ~ weibull (alpha=1/sigma, lambda = e^-{X*beta}), lamda is called scale
  #log(t_i*) ~ Gumble ( e^{X*beta}, sigma)
  #log(t_i*) = X*beta + sigma * epsilon, where epsilon ~ Gumble(0,1)
  sigma1 <- 1/2
  sigma0 <- 1/2
  alpha1 <- 1/sigma1
  alpha0 <- 1/sigma0
  lp_a1 <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
  lp_a0 <- eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
  lambda1 <- exp(-lp_a1)
  lambda0 <- exp(-lp_a0)
  D_a1 <- rweibull(n, shape = alpha1 , scale = 1/lambda1)
  D_a0 <- rweibull(n, shape = alpha0 , scale = 1/lambda0)
  #2nd method to generate D
  #Gumble(0,1)
  # epsilon1 <- -log(-log(runif(n)))  # for D_a1
  # epsilon0 <- -log(-log(runif(n)))  # for D_a0
  # 
  # D_a1 <- exp(lp_a1 + sigma * epsilon1)
  # D_a0 <- exp(lp_a0 + sigma * epsilon0)
  
  #CDF for weibull aft is: F(t)=1-exp(-lamda * t^sigma) -> T =  ((-log(-(1-F)))/lambda)^(1/sigma)
  # random_simu <- runif(n) #1-F(t) ~ Unif(0,1)
  # shape = 2 #weibull with shape = 1/sigma
  # lambda = 0.1 #base line rate, T~Weibull(v, rate) with rate for each individual is lambda*e^(X*beta)
  #T_i = (-log(U)/ lambda_i)^(1/v), lambda_i = lambda * e^(X*beta)
  #survival time D=(-log(random_simu) / (lambda * exp(eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U))) ^ (1 / v)
  # D_a0 = (-log(random_simu) / (lambda * exp(eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U))) ^ (1 / v)
  # D_a1 = (-log(random_simu) / (lambda * exp(eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U))) ^ (1 / v)

  
  # Observed time M=min(C,D,tau)
  M <- ifelse(A == 1, pmin(tau, C_a1, D_a1), pmin(tau, C_a0, D_a0))
  
  # Event indicator, delta=1 for event happened, otherwise delta=0
  delta <- ifelse(M == tau, 0,
                  ifelse(A == 1, as.numeric(D_a1 <= C_a1), as.numeric(D_a0 <= C_a0)))
  
  # Censoring indicator, gamma=1 for censored, otherwise gamma=0
  gamma <- ifelse(M == tau, 0,
                  ifelse(A == 1, as.numeric(C_a1 < D_a1), as.numeric(C_a0 < D_a0)))
  
  # simulated dataset (observed data)
  data_sim <- data.frame(
    M = M,
    status = ifelse(delta == 1, 1, 0),
    A = A,
    x1 = X1,
    x2 = X2,
    U = U
  )
  
  #compare
  data_sim_com <- data.frame(
    D_a0 = D_a0, D_a1=D_a1, C_a0 =C_a0, C_a1 = C_a1,
    M = M,
    status = ifelse(delta == 1, 1, 0),
    A = A,
    x1 = X1,
    x2 = X2,
    U = U
  )
  
  #Point estimate
  spce_withU[i] <- compute_spce_withU(data_sim)
  spce_withoutU[i] <- compute_spce_withoutU(data_sim) 
  
  #IPTW with/without U fpr SPCE
  #questions??? IPTW is nonparametric for marginal potential survival curves- KM estimator
  source("IPTW.R")
  #ipw_spce_withU(data, time_var, event_var, treat_var, t_eval)
  ipw_spce_withU[i] <- compute_ipw_spce_withU(data = data_sim,time_var = "M",event_var = "status",
                                              treat_var = "A",t_eval = t_pred)$SPCE
  ipw_spce_withoutU[i] <- compute_ipw_spce_withoutU(data = data_sim,time_var = "M",event_var = "status",
                                                    treat_var = "A",t_eval = t_pred)$SPCE
  
  
  
  #use Huang's function
  # Step 1: copy Huang's 2 functions for STOEM_IPW to working directory;
  source("SimulateU_surv_wb_2.R")
  source("surv_stoEM_ipw_wb_2.R")
  # # Step 2: update STOEM_IPW with saving probabiliy over HR, justing simple predict (HRmdoel, se.xxx=True);
  # # Apply your data_sim to their function;
  spce_huang[i] <- surv_stoEM_ipw(data=data_sim, zetat=-0.8, zetaz=-0.8, B = 50, theta = 0.5)$spce
  
  
  #bootstrap Within-Simulation sd
  for(b in 1:nboot){
    boot.index <- sample(x=1:n, size = n, replace = T)
    data.boot <- data_sim[boot.index, ]
    
    spce_withU.b[b] <- compute_spce_withU(data.boot)
    spce_withoutU.b[b] <- compute_spce_withoutU(data.boot)
    ipw_spce_withU.b[b] <- compute_ipw_spce_withU(data = data.boot,time_var = "M",event_var = "status",
                                                  treat_var = "A",t_eval = t_pred)$SPCE
    ipw_spce_withoutU.b[b] <- compute_ipw_spce_withoutU(data = data.boot,time_var = "M",event_var = "status",
                                                     treat_var = "A",t_eval = t_pred)$SPCE
    
    spce_huang.b[b] <- surv_stoEM_ipw(data=data_sim, zetat=-0.8, zetaz=-0.8, B = 50, theta = 0.5)$spce
    
  }
  
  
  spce_withU.sd[i] <- sd(spce_withU.b)
  spce_withoutU.sd[i] <-sd(spce_withoutU.b)
  spce_withU_iptw.sd[i] <- sd(ipw_spce_withU.b)
  spce_withoutU_iptw.sd[i] <- sd(ipw_spce_withoutU.b)
  spce_huang.sd[i] <- sd(spce_huang.b)
  
  
  cat("Simulation", i, "is finished\n")
}

#cp 
#G-computation with U
ci.lower.withU<-  spce_withU - qnorm(0.975) * spce_withU.sd
ci.upper.withU <-  spce_withU + qnorm(0.975) * spce_withU.sd
cp_withU <- true_spce >= ci.lower.withU & true_spce <= ci.upper.withU
mean(cp_withU)
#G-computation without U
ci.lower.withoutU <-  spce_withoutU - qnorm(0.975) * spce_withoutU.sd
ci.upper.withoutU <-  spce_withoutU + qnorm(0.975) * spce_withoutU.sd
cp_withoutU <- true_spce >= ci.lower.withoutU & true_spce <= ci.upper.withoutU
mean(cp_withoutU)

#iptw with U
ci.lower.iptw.withU <-  ipw_spce_withU - qnorm(0.975) * spce_withU_iptw.sd
ci.upper.iptw.withU <-  ipw_spce_withU + qnorm(0.975) * spce_withU_iptw.sd
cp_iptw.withU <- true_spce >= ci.lower.iptw.withU & true_spce <= ci.upper.iptw.withU
mean(cp_iptw.withU)
#iptw without U
ci.lower.iptw.withoutU <-  ipw_spce_withoutU - qnorm(0.975) * spce_withoutU_iptw.sd
ci.upper.iptw.withoutU <-  ipw_spce_withoutU + qnorm(0.975) * spce_withoutU_iptw.sd
cp_iptw.withoutU <- true_spce >= ci.lower.iptw.withoutU & true_spce <= ci.upper.iptw.withoutU
mean(cp_iptw.withoutU)

ci.lower.huang <-  spce_huang - qnorm(0.975) * spce_huang.sd
ci.upper.huang <-  spce_huang + qnorm(0.975) * spce_huang.sd
cp_huang <- true_spce >= ci.lower.huang & true_spce <= ci.upper.huang
mean(cp_huang)

#mean
mean(spce_withU[i])
mean(spce_withoutU[i])
mean(ipw_spce_withU[i])
mean(ipw_spce_withoutU[i])
mean(spce_huang[i])