#Simulation
#G-computation method
library(survival)
#set.seed(2025)
n_sim <-3 #number of simulation
n <- 1000 #sample size 
tau <- 5.5 #maximum follow-up time
t_pred <- 2  #time point for estimates 
true_spce <- 0.1482882#true SPCE at time=2 from numeric iteration
nboot<-100
spce_withU <- numeric(n_sim)
spce_withoutU <- numeric(n_sim)
spce_huang <- numeric(n_sim)
spce_withU.b <- rep(NA, nboot)
spce_withoutU.b <- rep(NA, nboot)
spce_huang.b <- rep(NA, nboot)
spce_withU.sd <- numeric(n_sim)
spce_withoutU.sd <- numeric(n_sim)
spce_huang.sd <- numeric(n_sim)
cp_withU <- logical(n_sim)
cp_withoutU <- logical(n_sim)
cp_huang <- logical(n_sim)


#function for computing SPCE
compute_spce_withU <- function(data) {
  # Fit seperated log-normal AFT models for A=1 and A=0 by using sub-dataset
  model1 <- coxph(Surv(M, status) ~ x1 + x2 + U, data = subset(data, A == 1))
  model0 <- coxph(Surv(M, status) ~ x1 + x2 + U, data = subset(data, A == 0))
  # Predict survival probability at specific time point over entire dataset
  mu1 <- predict(model1, newdata = data, type = "lp")
  mu0 <- predict(model0, newdata = data, type = "lp")
  
  S1 <- predict(model1, newdata=data, type = "survival", times = t_pred)
  S0 <- predict(model0, newdata=data, type = "survival", times = t_pred)
  # S1 <- 1 - pnorm((log(t_pred) - mu1) / model1$scale, mean=0, sd=1) 
  # S0 <- 1 - pnorm((log(t_pred) - mu0) / model0$scale, mean=0, sd=1) 
  #S(t)= 1- phi((log(t)-mu)/sigma), phi is CDF of standard normal
  #OR S <- 1 - plnorm(t_pred, meanlog = mu, sdlog = model$scale) 
  
  return(mean(S1) - mean(S0))# Estimate SPCE
}

compute_spce_withoutU <- function(data) {
  # Fit seperated log-normal AFT models for A=1 and A=0 by using sub-dataset
  # model1 <- survreg(Surv(M, status) ~ x1 + x2 , data = subset(data, A == 1), dist = "lognormal")
  # model0 <- survreg(Surv(M, status) ~ x1 + x2 , data = subset(data, A == 0), dist = "lognormal")
  # # Predict survival probability at specific time point over entire dataset
  # mu1 <- predict(model1, newdata = data, type = "lp")
  # mu0 <- predict(model0, newdata = data, type = "lp")
  # 
  # S1 <- 1 - pnorm((log(t_pred) - mu1) / model1$scale, mean=0, sd=1) 
  # S0 <- 1 - pnorm((log(t_pred) - mu0) / model0$scale, mean=0, sd=1) 
  # #S(t)= 1- phi((log(t)-mu)/sigma), phi is CDF of standard normal
  # #OR S <- 1 - plnorm(t_pred, meanlog = mu, sdlog = model$scale) 
  
  model1 <- coxph(Surv(M, status) ~ x1 + x2, data = subset(data, A == 1))
  model0 <- coxph(Surv(M, status) ~ x1 + x2, data = subset(data, A == 0))
  # Predict survival probability at specific time point over entire dataset
  mu1 <- predict(model1, newdata = data, type = "lp")
  mu0 <- predict(model0, newdata = data, type = "lp")
  
  S1 <- predict(model1, newdata=data, type = "survival", times = t_pred)
  S0 <- predict(model0, newdata=data, type = "survival", times = t_pred)
  
  return(mean(S1) - mean(S0)) # Estimate SPCE
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
  # A=0 parameters
  #censored
  # delta_intercept_a0 <- 0.6
  # delta_x1_a0 <- 0.3
  # delta_x2_a0 <- 0.7
  # delta_u_a0 <- -0.1
  # #event
  # eta_intercept_a0 <- 0.5
  # eta_x1_a0 <- -0.3
  # eta_x2_a0 <- 0.2
  # eta_u_a0 <- -1
  # # A=1 parameters
  # #censored
  # delta_intercept_a1 <- 1
  # delta_x1_a1 <- 0.7
  # delta_x2_a1 <- 1
  # delta_u_a1 <- 0.2
  #event
  eta_intercept <- 0.7
  eta_x1 <- -0.1
  eta_x2 <- 0.4
  eta_u <- -0.8
  eta_intercept_a1<- 0.2
  # 
  #   #C_a0 potential censored time for A=0
  #   C_a0 <- exp(delta_intercept_a0 + delta_x1_a0 * X1 + delta_x2_a0 * X2 + delta_u_a0 * U + rnorm(n, 0, 0.4))
  #   #D_a0 potential event time for A=0
  #   D_a0 <- exp(eta_intercept_a0 + eta_x1_a0 * X1 + eta_x2_a0 * X2 + eta_u_a0 * U + rnorm(n, 0, 0.4))
  #   #C_a1 potential censored time for A=1
  #   C_a1 <- exp(delta_intercept_a1 + delta_x1_a1 * X1 + delta_x2_a1 * X2 + delta_u_a1 * U + rnorm(n, 0, 0.4))
  #   #D_a1 potential event time for A=1
  #   D_a1 <- exp(eta_intercept_a1 + eta_x1_a1 * X1 + eta_x2_a1 * X2 + eta_u_a1 * U + rnorm(n, 0, 0.4))
  
  
  #C_a0 potential censored time for A=0
  C_a0 <- C_a1 <- runif(n, 0.1, 5.5)
  #D_a0 potential event time for A=0
  # D_a0 <- exp(eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U + rnorm(n, 0, 0.4))
  #C_a1 potential censored time for A=1
  #D_a1 potential event time for A=1
  # D_a1 <- exp(eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U + rnorm(n, 0, 0.4))
  
  random_simu <- runif(n)
  v =2
  lambda = 0.1
  D_a0 = round((-log(random_simu) / (lambda * exp(eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U))) ^ (1 / v))
  D_a1 = round((-log(random_simu) / (lambda * exp(eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U))) ^ (1 / v))
  
  
  # Survival_time = Survival_time1*Z + Survival_time0*(1-Z)
  
  
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
  
  #Point estimate
  
  spce_withU[i] <- compute_spce_withU(data_sim)
  spce_withoutU[i] <- compute_spce_withoutU(data_sim)
  
  #use Huang's function
  # Step 1: copy Huang's 2 functions for STOEM_IPW to working directory;
  source("SimulateU_surv.R")
  source("surv_stoEM_ipw.R")
  # Step 2: update STOEM_IPW with saving probabiliy over HR, justing simple predict (HRmdoel, se.xxx=True);
  # Apply your data_sim to their function;
  spce_huang[i] <- surv_stoEM_ipw(data=data_sim, zetat=-0.8, zetaz=-0.8, B = 50, theta = 0.5)$spce
  
  
  # #bootstrap
  # for(b in 1:nboot){
  #   boot.index <- sample(x=1:n, size = n, replace = T)
  #   data.boot <- data_sim[boot.index, ]
  #   
  #   spce_withU.b[b] <- compute_spce_withU(data.boot)
  #   spce_withoutU.b[b] <- compute_spce_withoutU(data.boot)
  #   spce_huang.b[b] <- surv_stoEM_ipw(data=data_sim, zetat=-0.8, zetaz=-0.8, B = 10, theta = 0.5)$spce
  #   
  # }
  # 
  # 
  # 
  # #standard deviation
  # spce_withU.sd[i] <- sd(spce_withU.b)
  # spce_withoutU.sd[i] <-sd(spce_withoutU.b)
  # spce_huang.sd[i] <- sd(spce_huang.b)
  
  cat("Simulation", i, "is finished\n")
}

#cp 
ci.lower.withU<-  spce_withU - qnorm(0.975) * spce_withU.sd
ci.upper.withU <-  spce_withU + qnorm(0.975) * spce_withU.sd
cp_withU <- true_spce >= ci.lower.withU & true_spce <= ci.upper.withU
mean(cp_withU)

ci.lower.withoutU <-  spce_withoutU - qnorm(0.975) * spce_withoutU.sd
ci.upper.withoutU <-  spce_withoutU + qnorm(0.975) * spce_withoutU.sd
cp_withoutU <- true_spce >= ci.lower.withoutU & true_spce <= ci.lower.withoutU
mean(cp_withoutU)

ci.lower.huang <-  spce_huang - qnorm(0.975) * spce_huang.sd
ci.upper.huang <-  spce_huang + qnorm(0.975) * spce_huang.sd
cp_huang <- true_spce >= ci.lower.huang & true_spce <= ci.upper.huang
mean(cp_huang)

#mean
mean(spce_withU[i])
mean(spce_withoutU[i])
mean(spce_huang[i])
