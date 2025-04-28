#Numeric Iteration for computing the true causal effect
#pre-specify the true values of parameters from data generating process

eta_intercept <- 0.7
eta_intercept_a1<- 0.2
eta_x1 <- -0.1
eta_x2 <- 0.4
eta_u <- -0.8


#
v = 2
lambda = 0.1 # lambda = (lambda*^v) = 0.1

#True SPCE

# Bernoulli probabilities
p_x1 <- function(x1) ifelse(x1 == 1, 0.4, 0.6) #X1~Bern(0.4)
p_x2 <- function(x2) ifelse(x2 == 1, 0.6, 0.4) #X1~Bern(0.6)
p_u <- function(u) 0.5  # U ~ Bern(0.5)

# Survival functions #Weibull/cox model S(t)= exp(-(lambda* * t)^v *exp(lp))
#Let lambda*^v = lambda, then S(t)= exp(-lambda * (t^v) *exp(lp))
#lp for A=1, eta_intercept_a1 + eta_x1 * x1 + eta_x2 * x2 + eta_u * u
#lp for A=0, eta_intercept + eta_x1 * x1 + eta_x2 * x2 + eta_u * u
S <- function(t, x1, x2, u, treat) {
  eta <- if (treat == 1) eta_intercept_a1 else eta_intercept
  lp <- eta + eta_x1 * x1 + eta_x2 * x2 + eta_u * u
  return(exp(-lambda * (t^v) * exp(lp)))
}

# E_X [ S(t | X, U, A) ]
E_X_given_U <- function(t, u, treat) {
  sum_val <- 0
  for (x1 in c(0,1)) {
    for (x2 in c(0,1)) {
      prob_x <- p_x1(x1) * p_x2(x2)
      sum_val <- sum_val + S(t, x1, x2, u, treat) * prob_x
    }
  }
  return(sum_val)
}

# E_U [ E_X [ S(t | X, U, A) ] ]
marginal_survival <- function(t, treat) {
  E_X_given_U(t, u=0, treat) * p_u(0) +
    E_X_given_U(t, u=1, treat) * p_u(1)
}

# SPCE
true_SPCE <- function(t) {
  marginal_survival(t, 1) - marginal_survival(t, 0)
}

t_val <-2
# Output
spce_val <- true_SPCE(t_val)
cat("True SPCE at time =", t_val, "is", round(spce_val, 6), "\n")
