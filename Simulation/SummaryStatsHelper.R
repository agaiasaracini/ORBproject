

############## Functions to extract statistics from simulations ################

# Function to calculate Bias and MCSE
calculate_bias_mcse <- function(theta_hat, theta) {
  N_sim <- length(theta_hat)
  bias <- mean(theta_hat - theta)
  mcse <- sqrt(sum((theta_hat - mean(theta_hat))^2) / (N_sim * (N_sim - 1)))
  return(list(bias = bias, mcse = mcse))
}

# Function to calculate Empirical SE and MCSE
calculate_empirical_se_mcse <- function(theta_hat) {
  N_sim <- length(theta_hat)
  empirical_se <- sqrt(sum((theta_hat - mean(theta_hat))^2) / (N_sim - 1))
  mcse <- empirical_se / sqrt(2 * (N_sim - 1))
  return(list(empirical_se = empirical_se, mcse = mcse))
}

# Function to calculate MSE and MCSE
calculate_mse_mcse <- function(theta_hat, theta) {
  N_sim <- length(theta_hat)
  mse <- mean((theta_hat - theta)^2)
  mse_hat <- mean((theta_hat - mean(theta_hat))^2)
  mcse <- sqrt(sum((theta_hat - theta)^2 - mse_hat) / (N_sim * (N_sim - 1)))
  return(list(mse = mse, mcse = mcse))
}

# Function to calculate Coverage and MCSE
calculate_coverage_mcse <- function(theta_hat_low, theta_hat_up, theta) {
  # Remove NA values
  valid_indices <- !is.na(theta_hat_low) & !is.na(theta_hat_up)
  theta_hat_low <- theta_hat_low[valid_indices]
  theta_hat_up <- theta_hat_up[valid_indices]
  
  N_sim <- length(theta_hat_low)
  coverage <- mean(theta_hat_low <= theta & theta <= theta_hat_up)
  mcse <- sqrt(coverage * (1 - coverage) / N_sim)
  return(list(coverage = coverage, mcse = mcse))
}

# Function to calculate Power and MCSE
calculate_power_mcse <- function(theta_hat_low, theta_hat_up) {
  # Remove NA values
  init <- length(theta_hat_low)
  valid_indices <- !is.na(theta_hat_low) & !is.na(theta_hat_up)
  theta_hat_low <- theta_hat_low[valid_indices]
  theta_hat_up <- theta_hat_up[valid_indices]
  
  N_sim <- length(theta_hat_low)
  power <- mean(theta_hat_low <= 0 & theta_hat_up >= 0)
  mcse <- sqrt(power * (1 - power) / N_sim)
  return(list(power = power, mcse = mcse, valid = length(valid_indices)/init))
}


############# Summary statistics function ######################################

# Function to calculate summary statistics for a given method
calculate_summary <- function(theta_hats, theta_hats_low, theta_hats_up, theta) {
  
  #N_sim <- length(theta_hats)
  
  # Calculate bias and MCSE
  bias_mcse <- calculate_bias_mcse(theta_hats, theta)
  
  # Calculate empirical SE and MCSE
  emp_se_mcse <- calculate_empirical_se_mcse(theta_hats)
  
  # Calculate MSE and MCSE
  mse_mcse <- calculate_mse_mcse(theta_hats, theta)
  
  # Calculate Coverage and MCSE
  cov_mcse <- calculate_coverage_mcse(theta_hats_low, theta_hats_up, theta)
  
  # Calculate Power and MCSE
  power_mcse <- calculate_power_mcse(theta_hats_low, theta_hats_up)
  
  # Combine all the summary statistics into a list
  summary <- list(
    bias = bias_mcse$bias,
    bias_mcse = bias_mcse$mcse,
    emp_se = emp_se_mcse$empirical_se,
    emp_se_mcse = emp_se_mcse$mcse,
    mse = mse_mcse$mse,
    mse_mcse = mse_mcse$mcse,
    cov = cov_mcse$coverage,
    cov_mcse = cov_mcse$mcse,
    power = power_mcse$power,
    power_mcse = power_mcse$mcse,
    valid = power_mcse$valid
    
  )
  
  return(summary)
}



