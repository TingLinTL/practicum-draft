#This dataset is generated event time only, no censoring
set.seed(2025)
n <- 1000 #sample size 
tau <- 5.5 #maximum follow-up time

#Generate X1, X2, U
X1 <- rbinom(n, size = 1, prob = 0.4) #X1~Bern(0.4)
X2 <- rbinom(n, size = 1, prob = 0.6) #X2~Bern(0.6)
U <- rbinom(n, size = 1, prob = 0.5) #U~Bern(0.5), independent with X1 and X2

# Treatment assignment, logit(P(A=1|X1,X2,U))
alpha_0 <- 0.1; alpha_x1 <- 0.3; alpha_x2 <- 0.5; alpha_u <- -0.8
pi_A <- 1 / (1 + exp(-(alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U)))
A <- rbinom(n, size = 1, prob = pi_A)

#event
eta_intercept <- 0.2
eta_intercept_a1<- 0.7
eta_x1 <- -0.1
eta_x2 <- 0.4
eta_u <- -0.8

#no censoring
#General weibull distribution
#log(T)=lp+sigma*ε, ε~Gumbel(0,1)
#T~Weibull(lamda(rate)=exp(-lp/sigma), v(shape)=1/sigma)
sigma <- 0.5 
mu0 <- eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U # Linear predictor 
mu1 <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
#transformation: log(T)=mu+sigma*error
# epsilon_star <- log(epsilon) ε*~ Gumbel(0,1),-log(-log(runif(n)))
u <- runif(n)          # Uniform(0,1)
gumbel_errors <- log(-log(u))  #Transform
# Generate survival time from Weibull AFT model
D_a0 <- exp(mu0+sigma*gumbel_errors)
D_a1 <- exp(mu1+sigma*gumbel_errors)

# Observed time M=min(C,D,tau)
M <- ifelse(A == 1, D_a1, D_a0)

#status, all time points are event time, which means status = 1
status<-rep(1,1000)

# simulated dataset (observed data)
data_sim <- data.frame(
  M = M,
  A = A,
  status = status,
  x1 = X1,
  x2 = X2,
  U = U
)
