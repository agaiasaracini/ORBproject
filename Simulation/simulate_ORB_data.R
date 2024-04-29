
############## Simulate meta-analysis dataset with ORB #########################

simulate_ORB_data <- function(n_studies,
                              mu,
                              tau_squared,
                              n_treatment=50,
                              n_control=50,
                              gamma=1.5) {
  
  # Function to simulate meta-analysis dataset with ORB
  simulate_single_dataset <- function(n_studies, mu, tau_squared, n_treatment, n_control, gamma) {
    # Empty data frame
    meta_data <- data.frame(study = 1:n_studies,
                            n_control = numeric(n_studies),
                            n_treatment = numeric(n_studies),
                            y = numeric(n_studies),
                            s = numeric(n_studies),
                            y_ORB = numeric(n_studies),
                            s_ORB = numeric(n_studies),
                            p_value = numeric(n_studies))
    
    # Generate true treatment effect mean and standard errors for all studies at once
    theta_i <- rnorm(n_studies, mean = mu, sd = sqrt(tau_squared))
    sigma_i <- sqrt(rchisq(n_studies, df = 2 * n_control - 2) / ((n_control - 1) * n_control))
    
    # Generate observed treatment effect for all studies at once
    y_i <- rnorm(n_studies, mean = theta_i, sd = sqrt(2 / n_control))
    
    # Calculate one-sided p-values for all studies at once
    p_i <- pnorm(y_i / sigma_i, lower.tail = FALSE)
    p_i <- round(p_i, 10)
    
    # Fill in data frame
    meta_data$n_control <- n_control
    meta_data$n_treatment <- n_treatment
    meta_data$y <- y_i
    meta_data$s <- sigma_i
    meta_data$y_ORB <- y_i
    meta_data$s_ORB <- sigma_i
    meta_data$p_value <- p_i
    
    
    
    # Selection process to simulate ORB
    replace_prob <- exp(-4 * as.numeric(meta_data$p_value)^gamma)
    unrep_index <- runif(n_studies) >= replace_prob
    meta_data$y_ORB[unrep_index] <- "unrep"
    meta_data$s_ORB[unrep_index] <- "unrep"
    
    return(meta_data)
  }
  
  # Initialize trial counter
  trials <- 0
  
  # Repeat the simulation until at least one study outcome is not "unrep"
  repeat {
    # Increment trial counter
    trials <- trials + 1
    
    # Generate meta-analysis dataset
    data <- simulate_single_dataset(n_studies, mu, tau_squared, n_treatment, n_control, gamma)
    
    # Check if there's at least one study outcome that is not "unrep"
    if (!all(data$y_ORB == "unrep")) {
      break  # Exit the loop if condition met
    }
  }
  
  num_rep <- n_studies - sum(data$y_ORB == "unrep")
  
  
  return(list(dataORB = data,
              Ntrials = trials,
              Nrep = num_rep))
}


