# Parallelized simulation using doParallel + foreach

#–– Enable JIT compilation for speed ––#
compiler::enableJIT(3)

#–– Load libraries ––#
library(parallel)
library(survival)

#–– Simulation parameters ––#
n_sim       <- 100          # number of simulati
n           <- 300          # sample size
tau         <- 5.5          # maximum follow-up time
t_pred      <- 2            # time point for SPCE estimation
true_spce   <- 0.1482882     # true SPCE at t_pred
dist_lambda <- 0.1          # Weibull scale
dist_v      <- 2            # Weibull shape
nboot       <- 100          # bootstrap replicates

#–– Define SPCE computation functions and byte-compile them ––#
compute_spce_withU <- compiler::cmpfun(function(data) {
  m1 <- coxph(Surv(M, status) ~ x1 + x2 + U, data = subset(data, A == 1))
  m0 <- coxph(Surv(M, status) ~ x1 + x2 + U, data = subset(data, A == 0))
  S1 <- predict(m1, newdata = data, type = "survival", times = t_pred)
  S0 <- predict(m0, newdata = data, type = "survival", times = t_pred)
  mean(S1) - mean(S0)
})
compute_spce_withoutU <- compiler::cmpfun(function(data) {
  m1 <- coxph(Surv(M, status) ~ x1 + x2, data = subset(data, A == 1))
  m0 <- coxph(Surv(M, status) ~ x1 + x2, data = subset(data, A == 0))
  S1 <- predict(m1, newdata = data, type = "survival", times = t_pred)
  S0 <- predict(m0, newdata = data, type = "survival", times = t_pred)
  mean(S1) - mean(S0)
})

#–– Define one-iteration simulation function and byte-compile ––#
sim_one <- compiler::cmpfun(function(i) {
  # 1. Covariates
  X1 <- rbinom(n, 1, 0.4)
  X2 <- rbinom(n, 1, 0.6)
  U  <- rbinom(n, 1, 0.5)
  
  # 2. Treatment
  lp <- 0.1 + 0.3*X1 + 0.5*X2 - 0.8*U
  A  <- rbinom(n, 1, plogis(lp))
  
  # 3. Event & censoring
  u   <- runif(n)
  D0  <- round((-log(u) / (dist_lambda * exp(0.7 - 0.1*X1 + 0.4*X2 - 0.8*U)))^(1/dist_v))
  D1  <- round((-log(u) / (dist_lambda * exp(0.2 - 0.1*X1 + 0.4*X2 - 0.8*U)))^(1/dist_v))
  C   <- runif(n, 0.1, 5.5)
  
  # 4. Observed
  M      <- ifelse(A == 1, pmin(tau, C, D1), pmin(tau, C, D0))
  status <- ifelse(M == tau, 0,
                   ifelse(A == 1, as.numeric(D1 <= C), as.numeric(D0 <= C)))
  
  dat <- data.frame(M = M, status = status, A = A, x1 = X1, x2 = X2, U = U)
  
  # 5. Estimates
  e1 <- compute_spce_withU(dat)
  e2 <- compute_spce_withoutU(dat)
  e3 <- surv_stoEM_ipw(dat, zetat = -0.8, zetaz = -0.8, B = 20, theta = 0.5)$spce
  
  # 6. Bootstrap SDs
  b1 <- numeric(nboot); b2 <- numeric(nboot); b3 <- numeric(nboot)
  for(b in seq_len(nboot)) {
    samp  <- sample.int(n, n, replace = TRUE)
    dboot <- dat[samp, ]
    b1[b] <- compute_spce_withU(dboot)
    b2[b] <- compute_spce_withoutU(dboot)
    b3[b] <- surv_stoEM_ipw(dboot, zetat = -0.8, zetaz = -0.8, B = 20, theta = 0.5)$spce
  }
  sd1 <- sd(b1); sd2 <- sd(b2); sd3 <- sd(b3)
  
  # 7. Coverage
  cp1 <- abs(e1 - true_spce) <= 1.96 * sd1
  cp2 <- abs(e2 - true_spce) <= 1.96 * sd2
  cp3 <- abs(e3 - true_spce) <= 1.96 * sd3
  
  c(e1, e2, e3, sd1, sd2, sd3, cp1, cp2, cp3)
})

#–– Parallel execution with load balancing ––#
cores <- detectCores(logical = FALSE) - 1
cl    <- makeCluster(cores)

# Export data and functions
clusterExport(cl, varlist = c(
  "n","tau","t_pred","true_spce",
  "dist_lambda","dist_v","nboot",
  "compute_spce_withU","compute_spce_withoutU",
  "sim_one"
))
# Ensure survival and IPW functions are available
clusterEvalQ(cl, library(survival))
clusterEvalQ(cl, {
  source("SimulateU_surv.R")    # defines SimulateU_surv and any dependencies
  source("surv_stoEM_ipw.R")    # defines surv_stoEM_ipw
})

# Run with load balancing
res_list <- parLapplyLB(cl, 1:n_sim, sim_one)
stopCluster(cl)

# Collect results
results_df <- do.call(rbind, res_list)
colnames(results_df) <- c(
  "spce_withU","spce_withoutU","spce_huang",
  "spce_withU_sd","spce_withoutU_sd","spce_huang_sd",
  "cp_withU","cp_withoutU","cp_huang"
)