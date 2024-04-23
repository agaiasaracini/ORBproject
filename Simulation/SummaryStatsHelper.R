

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
calculate_coverage_mcse <- function(theta_low, theta_up, theta) {
  N_sim <- length(theta_low)
  coverage <- mean(theta_low <= theta & theta <= theta_up)
  mcse <- sqrt(coverage * (1 - coverage) / N_sim)
  return(list(coverage = coverage, mcse = mcse))
}

# Function to calculate Power and MCSE
calculate_power_mcse <- function(p_values) {
  N_sim <- length(p_values)
  power <- mean(p_values <= 0.05)
  mcse <- sqrt(power * (1 - power) / N_sim)
  return(list(power = power, mcse = mcse))
}
