
# ------------------------------------------------------------------------------
# Data generating process and candidate model to estimate
# ------------------------------------------------------------------------------

# Data generation process
DGP <- "A"

# Candidate model
# (can be found in "~/Meta_DLC_SEM/candidate_models/model_specification/AR1)
candidate_model <- "AR1_no_crossloadings.txt"

# Random seed
random_seed <- 24

# ------------------------------------------------------------------------------
# Dimensions for the generated datasets and candidate models
# ------------------------------------------------------------------------------

# Number of persons
n_persons <- 100

# Number of time points
n_timepoints <- 25 

# Number of latent factors
n_factors <- 2

# Number of indicators for the latent factors
n_indicators <- 6

# Number of states the intra-individual processes can be in
n_states <- 2

# ------------------------------------------------------------------------------
# MCMC sampling parameters
# ------------------------------------------------------------------------------

# Number of iterations
n_iterations <- 10000

# Number of iterations used as warm-up period
n_warmup <- n_iterations / 2

# Number of chains
n_chains <- 5

# Mean and precision priors for the between-level latent variables
mu_prior <- c(0, 0)
precision_prior <- diag(2)

# Parameters estimates to retrieve
MCMC_params <- c(
  
  # AR(1) intercepts and coefficients
  "alpha_S1", "alpha_S2", "beta_within_S1", "beta_within_S2",
  
  # Factor loadings
  "lambda_within_S1", "lambda_within_S2",
  
  # State memberships and time point of state switch
  "state_membership", "state_switch",
  
  # Predicted indicator scores
  "mu_y_within"
  
)

# ------------------------------------------------------------------------------
# Fixed parameter values for the generated datasets
# ------------------------------------------------------------------------------

# The parameters are defined for the case of two latent factors and states;
# for more factors, the variables need to be adjusted accordingly

# Mean vector for intra-individual differences
mu_RE <- c(0, 0)

# Factor variances
factor_variance <- 0.3

# Residual variances on the between- and within-level are identical
residual_variance <- 0.25

# State-specific factor means
mu_S1 = c(0, 0)
mu_S2 = c(2, 2)

# State-specific AR(1) coefficients (intercepts are fixed at 0)
beta_S1 <- c(.6, .6)
beta_S2 <- c(.3, .3)

# State-probability
pi = 0.5

# Support
support <- matrix(1, n_persons, n_factors)
