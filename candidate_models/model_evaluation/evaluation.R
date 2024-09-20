
suppressPackageStartupMessages(
  
  {
    
    library(boot)
    library(glue)
    library(LaplacesDemon)
    library(loo)
    library(R2jags)
    
  }
  
)

options(scipen = 999)
set.seed(24)

# Generating the targets
source(glue("data_generation/generate_dataset_{DGP}.R"))

# Calibrating the candidate specifications in JAGS -----------------------------

# Path of the candidate model
model_file <- glue("candidate_models/model_specification/{candidate_model}.txt")

# Preparing the JAGS model
MCMC_data <- list(
  
  # Dimensions
  n_persons = n_persons,
  n_timepoints = n_timepoints,
  n_indicators = n_indicators,
  n_factors = n_factors, 
  n_states = n_states,
  
  # Taking the generated ground truth data as targets
  y = y,
  
  # Mean and precision prior
  mu_0 = mu_prior,
  phi_0 = precision_prior
  
)


# Calibrating the model
MCMC_results <- jags.parallel(
  
  data = MCMC_data, 
  parameters.to.save = MCMC_params,
  model.file = model_file,
  n.chains = n_chains,
  n.iter = n_iterations,
  n.burnin = n_warmup,
  export = c("n_chains", "n_iterations", "n_warmup") 
  
)

# Extracting model outputs -----------------------------------------------------

# Note that these might need to be adapted depending on the MCMC_params defined in configuration.R

# Score estimates
y_hat_raw <- MCMC_results$BUGSoutput$summary[
  (
    1 +
      1 +
      2 * (n_factors * n_states) + 
      2 * (n_indicators - n_factors)
  ):
    (
      0 + 
        1 + 
        2 * (n_factors * n_states) + 
        2 * (n_indicators - n_factors) +
        n_timepoints * n_persons * n_indicators * n_states
    ), 
  "50%"
]

# State membership estimates
state_membership_raw <- MCMC_results$BUGSoutput$summary[
  (
    1 +
      1 +
      2 * (n_factors * n_states) +
      2 * (n_indicators - n_factors) + 
      n_timepoints * n_persons * n_indicators * n_states 
  ):
    (
      0 + 
        1 +
        2 * (n_factors * n_states) + 
        2 * (n_indicators - n_factors) + 
        n_timepoints * n_persons * n_indicators * n_states + 
        n_persons * n_timepoints
    ), 
  "50%"
]

# State switch estimates
state_switch_hat <- MCMC_results$BUGSoutput$summary[
  (
    1 +                                                                         
      1 +                                                                       
      2 * (n_factors * n_states) +                                              
      2 * (n_indicators - n_factors) +                                          
      n_timepoints * n_persons * n_indicators * n_states +                      
      n_persons * n_timepoints                                                  
  ):
    (
      0 +
        1 +
        2 * (n_factors * n_states) +
        2 * (n_indicators - n_factors) + 
        n_timepoints * n_persons * n_indicators * n_states +
        n_persons * n_timepoints +
        n_persons
    ),
  "50%"
]

# Store state-specific score estimates
y_hat_S1 <- y_hat_raw[1:(length(y_hat_raw) / 2)] 
y_hat_S2 <- y_hat_raw[(1 + length(y_hat_raw) / 2):length(y_hat_raw)]

# Create a meaningful score estimate for loss calculation
y_hat <- array(NA, dim = c(n_timepoints, n_persons, n_indicators))
state_membership_matrix <- matrix(
  
  state_membership_raw, 
  nrow = n_timepoints, 
  ncol = n_persons, 
  byrow = T
  
)

for (i in 1:n_persons) {
  
  for (t in 1:n_timepoints) {
    
    state <- state_membership_matrix[t, i]
    idx_start <- ((i - 1) * n_timepoints * n_indicators) + 
      ((t - 1) * n_indicators) + 1
    idx_end <- idx_start + n_indicators - 1
    
    if (state == 1) {
      y_hat[t, i, ] <- y_hat_S1[idx_start:idx_end]
    } 
    
    else {
      y_hat[t, i, ] <- y_hat_S2[idx_start:idx_end]
    }
    
  }
  
}

# Storing calibration results --------------------------------------------------

# Calculating the loss and weight surrogates (multiple-steps ahead) 
(SE <- apply(((y_hat - y) ^ 2), c(1), sum))
(MSE <- SE / (n_persons * n_indicators))
(MSE_inv <- 1 / MSE)

# Storing results
write.csv(
  y_hat, 
  glue("candidate_models/model_output/predictions/pred_{candidate_model}_{DGP}.csv")
)
write.csv(
  MSE_inv, 
  glue("candidate_models/model_output/loss/loss_{candidate_model}_{DGP}.csv")
)
write.csv(
  state_switch_hat, 
  glue("candidate_models/model_output/state_switch/state_switch_{candidate_model}_{DGP}.csv")
)
