
############## Simulate meta-analysis dataset with ORB, MAR #########################

simulate_ORB_data_MAR <- function(n_studies,
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
    
    dataALL <- meta_data
    
    # Selection process to simulate ORB
    replace_prob <- exp(-4 * as.numeric(meta_data$p_value)^gamma)
    unrep_index <- runif(n_studies) >= replace_prob
    meta_data$y_ORB[unrep_index] <- "unrep"
    meta_data$s_ORB[unrep_index] <- "unrep"
    
    return(list(meta_data,
                dataALL))
  }
  
  # Initialize trial counter
  trials <- 0
  
  # Repeat the simulation until at least one study outcome is not "unrep"
  repeat {
    # Increment trial counter
    trials <- trials + 1
    
    # Generate meta-analysis dataset
    d <- simulate_single_dataset(n_studies, mu, tau_squared, n_treatment, n_control, gamma)
    data_noORB <- d[[2]]
    data <- d[[1]]
    
    # Check if there's at least two studies outcomes that is not "unrep"
    if (sum(!(data$y_ORB == "unrep")) > 2) {
      break  # Exit the loop if condition met
    }
  }
  
  num_rep <- n_studies - sum(data$y_ORB == "unrep")
  
  num_unrep <- sum(data$y_ORB == "unrep")
  
  # Create dataMAR: same data, but remove 'num_unrep' at random
  dataMAR <- data_noORB
  rand_indices <- sample(1:n_studies, num_unrep, replace = FALSE)
  dataMAR$y_ORB[rand_indices] <- "unrep"
  dataMAR$s_ORB[rand_indices] <- "unrep"
  
  
  return(list(dataORB = dataMAR,
              Ntrials = trials,
              Nrep = num_rep))
}


