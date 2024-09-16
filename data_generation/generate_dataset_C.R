
library(mvtnorm)

# Data is generated based on previously defined dimensions and parameter values
source("~/Meta_DLC_SEM/configuration.R")

options(scipen = 999)
set.seed(random_seed)

# Factor loading matrices ------------------------------------------------------

# In state 1 and prior to time point 15, indicators 1, 2, and 3 load on factor 1,
# the remaining load on factor 2; 
# after 15 points in time indicator 4 additionally loads on factor 1;
# in state 2 indicators 1, 2, and 3 load on factor 1, the remaining load on factor 2

# State 1 for t = 1,  ..., 15
lambda_S1_before <- matrix(
  
  c(rep(1, (n_indicators / n_factors)), rep(0, (n_indicators / n_factors)), 
    rep(0, (n_indicators / n_factors)), rep(1, (n_indicators / n_factors))),
  n_indicators, n_factors
  
)

# State 1 for t = 16, ..., n_timepoints
lambda_S1_after <- matrix(
  
  c(rep(1, 4), rep(0, 2), 
    rep(0, (n_indicators / n_factors)), rep(1, (n_indicators / n_factors))),
  n_indicators, n_factors
  
)

# State 2 (factor loadings are time-invariant)
lambda_S2 <- matrix(
  
  c(rep(1, 4), rep(0, 2), 
    rep(0, (n_indicators / n_factors)), rep(1, (n_indicators / n_factors))),
  n_indicators, n_factors
  
)

# State-specific variances -----------------------------------------------------

# Residual variances
epsilon_S1 <- diag(n_indicators) * residual_variance
epsilon_S2 <- diag(n_indicators) * residual_variance

# Factor variances
phi_within_S1 <- diag(2)
phi_between_S1 <- phi_within_S1 * factor_variance
phi_within_S2 <- diag(2)
phi_between_S2 <- phi_within_S2 * factor_variance

# Factor and indicator scores --------------------------------------------------

# Between-level random effects are identical in both states
eta_between_S1 <- eta_between_S2 <- rmvnorm(n_persons, mu_RE, phi_between_S1)

# State 1:
# --------

# Initialize within-level factor scores
eta_within_S1 <- array(NA, c(n_timepoints, n_persons, n_factors))

# Initialize indicator scores
y_S1 <- array(NA, c(n_timepoints, n_persons, n_indicators))

# Determine scores through time

# t = 1
eta_within_S1[1, , ] <- mu_S1 * support + eta_between_S1 + rmvnorm(
  n_persons, mu_RE, phi_within_S1 - phi_between_S1
)
y_S1[1, , ] <- eta_within_S1[1, , ] %*% t(lambda_S1) + rmvnorm(
  n_persons, sigma = epsilon_S1
)

# t = 2, ..., 15
for (t in 2:15){
  
  # Residuals
  zeta_within <- rmvnorm(
    n_persons, mu_RE, phi_within_S1 * (1 - beta_S1 ^ 2) - phi_between_S1
  )
  
  # Factor scores
  eta_within_S1[t, , 1] <- mu_S1[1] + eta_between_S1[, 1] + beta_S1[1] * (
    eta_within_S1[t - 1, , 1] - eta_between_S1[, 1] - mu_S1[1]
  ) + zeta_within[, 1]
  eta_within_S1[t, , 2] <- mu_S1[2] + eta_between_S1[, 2] + beta_S1[2] * (
    eta_within_S1[t - 1, , 2] - eta_between_S1[, 2] - mu_S1[2]
  ) + zeta_within[, 2]
  
  # Observed scores
  y_S1[t, , ] <- eta_within_S1[t, , ] %*% t(lambda_S1_before) + rmvnorm(
    n_persons, sigma = epsilon_S1
  )
  
}

# t = 16, ..., n_timepoints
for (t in 16:n_timepoints){
  
  # Residuals
  zeta_within <- rmvnorm(
    n_persons, mu_RE, phi_within_S1 * (1 - beta_S1 ^ 2) - phi_between_S1
  )
  
  # Factor scores
  eta_within_S1[t, , 1] <- mu_S1[1] + eta_between_S1[, 1] + beta_S1[1] * (
    eta_within_S1[t - 1, , 1] - eta_between_S1[, 1] - mu_S1[1]
  ) + zeta_within[, 1]
  eta_within_S1[t, , 2] <- mu_S1[2] + eta_between_S1[, 2] + beta_S1[2] * (
    eta_within_S1[t - 1, , 2] - eta_between_S1[, 2] - mu_S1[2]
  ) + zeta_within[, 2]
  
  # Observed scores
  y_S1[t, , ] <- eta_within_S1[t, , ] %*% t(lambda_S1_after) + rmvnorm(
    n_persons, sigma = epsilon_S1
  )
  
}

# State 2:
# --------

# Initialize within-level factor scores
eta_within_S2 <- array(NA, c(n_timepoints, n_persons, n_factors))

# Initialize indicator scores
y_S2 <- array(NA, c(n_timepoints, n_persons, n_indicators))

# Determine scores through time

# t = 1
eta_within_S2[1, , ] <- mu_S2 * support + eta_between_S2 + rmvnorm(
  n_persons, mu_RE, phi_within_S2 - phi_between_S2
)
y_S2[1, , ] <- eta_within_S2[1, , ] %*% t(lambda_S2_before) + rmvnorm(
  n_persons, sigma = epsilon_S2
)

# t = 2, ..., n_timepoints
for (t in 2:n_timepoints){
  
  # Residuals
  zeta_within <- rmvnorm(
    n_persons, mu_RE, phi_within_S2 * (1 - beta_S2 ^ 2) - phi_between_S2
  )
  
  # Factor scores
  eta_within_S2[t, , 1] <- mu_S2[1] + eta_between_S2[, 1] + beta_S2[1] * (
    eta_within_S2[t - 1, , 1] - eta_between_S2[, 1] - mu_S2[1]
  ) + zeta_within[, 1]
  eta_within_S2[t, , 2] <- mu_S2[2] + eta_between_S2[, 2] + beta_S2[2] * (
    eta_within_S2[t - 1, , 2] - eta_between_S2[, 2] - mu_S2[2]
  ) + zeta_within[, 2]
  
  # Indicator scores
  y_S2[t, , ] <- eta_within_S2[t, , ] %*% t(lambda_S2) + rmvnorm(
    n_persons, sigma = epsilon_S2
  )
  
}

# State memberships ------------------------------------------------------------

# Define switching individuals by sampling with replacement and probability pi
state_membership <- sample(
  
  c(1, 2), 
  n_persons, 
  prob = c(pi, 1 - pi), 
  replace = T
  
)

# Each individual starts in state 1
y <- y_S1

# Switching individuals switch at their time of state switch and stay in state 2
state_switch <- sample(1:n_timepoints, n_persons, replace = T)
for (i in 1:n_persons){
  
  if(state_membership[i] == 2){
    
    y[state_switch[i]:n_timepoints, i, ] <- y_S2[
      state_switch[i]:n_timepoints, i,
    ]
    
  }
  
}

# Results are saved
write.csv(y, "~/Meta_DLC_SEM/data_generation/target/y_C.csv")
write.csv(state_switch, "~/Meta_DLC_SEM/data_generation/state_switch/state_switch_C.csv")
