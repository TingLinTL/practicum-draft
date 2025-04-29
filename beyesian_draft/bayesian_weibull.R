options(scipen = 999)
library(coda)
library(rjags)
library(dplyr)
library(survival)

source("Baye_weibull_data generating_noC.R")

cat("
model {
  for (i in 1:N) {
  
    # P(T,A,U|X)=P(T|A,X,U)*P(A|X,U)*P(U|X)
    
    #  Likelihood for U (latent confounder), P(U|X), if U is dependent of X1 and X2
    U[i] ~ dbern(pU[i])
    logit(pU[i])<- alpha0 + alpha.X1 * X1[i]+ alpha.X2 * X2[i]
    
    # Treatment A, P(A|X,U)
    A[i] ~ dbern(pA[i])
    logit(pA[i]) <- gamma0 + gamma.X1 * X1[i] + gamma.X2 * X2[i]+ gamma.U * U[i]
    
    # Likelihood for survival outcome T (AFT model), P(T|A,X,U)
    
    #the following contents are needed to be verified again
    #lambda is the rate parameter, not scale, so in dweib() of JAGS, lamda = exp(-mu[i]/sigma)
    # T~weibull(shape=1/sigma, rate= exp(-mu[i]/sigma))
    #T[i] ~ dweib(shape, lambda[i])   Weibull survival time, dweib(shape, rate)
    #lambda[i] <- exp(-mu[i]/sigma)
    
    
    #Weibull logT~N(mu, tau)
    logT[i] ~ dnorm(mu[i], tau) 
    mu[i] <- beta0 + beta.A * A[i] + beta.X1 * X1[i] + beta.X2 * X2[i] + beta.U * U[i]
    
    # spce;
    # U is independent of A
    # predicted potential mu under A=1 and A=0
    pUpo[i] <- ilogit(alpha0 + alpha.X1 * X1[i]+ alpha.X2 * X2[i])
    mu1po[i] <- beta0 + beta.A * 1 + beta.X1 * X1[i] + beta.X2 * X2[i] + beta.U * pUpo[i]
    mu0po[i] <- beta0 + beta.A * 0 + beta.X1 * X1[i] + beta.X2 * X2[i] + beta.U * pUpo[i]
    S1[i] <- 1 - phi((log(t_pred) - mu1po[i])/sigma)
    S0[i] <- 1 - phi((log(t_pred) - mu0po[i])/sigma)
    
  }
  
  spce <- mean(S1[])-mean(S0[])
  
  # Priors for parameters
  sigma ~ dunif(0.01, 10)
  tau <- 1/sigma/sigma
  
  gamma0 ~ dnorm(0, 0.01)
  gamma.X1 ~ dnorm(0, 0.01)
  gamma.X2 ~ dnorm(0, 0.01)
  
  beta0 ~ dnorm(0, 0.01)
  beta.A ~ dnorm(0, 0.01)
  beta.X1 ~ dnorm(0, 0.01)
  beta.X2 ~ dnorm(0, 0.01)
  
  #bias parameters
  gamma.U ~ dunif(-5,5)
  alpha0 ~ dunif(-5,5)
  alpha.X1 ~ dunif(-5,5)
  alpha.X2 ~ dunif(-5,5)
  beta.U ~ dunif(-5,5)

}
 ", file="weibull_JAGS.txt")


jags_data <- list(
  N = nrow(data_sim),
  logT = log(data_sim$M),
  X1 = data_sim$x1,
  X2 = data_sim$x2,
  A = data_sim$A,
  U = data_sim$U,
  t_pred = 2
)


# Model
model <- jags.model("weibull_JAGS.txt", data = jags_data, n.chains = 3, n.adapt = 5000)

update(model, 1000) # Burn-in

samples <- coda.samples(
  model, variable.names = c("beta0","beta.A", "beta.X1","beta.X2","beta.U", "sigma","spce"), 
  n.iter = 20000, thin = 10)


summary(samples)
mcmc <- as.matrix(samples)
traceplot(samples[, "spce"], main = "Traceplot of SPCE")
acf(mcmc[, "spce"], main = "ACF of SPCE")


#compare
source("spce functions.R")

t_pred <-2
compute_spce_withU(data_sim)
compute_spce_withoutU(data_sim)

