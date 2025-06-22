

#New code, all together for efficiency and memory usage

# Load necessary sources and libraries
source("simulate_ORB_data_MAR.R")
source("reORBadjust.R")
source("SummaryStatsHelper.R")

library(dplyr)
library(foreach)
library(doParallel)
library(doRNG)

# Set seed for reproducibility
seed <- 123

# Grid of parameters for the DGM
parDATA_MAR <- expand.grid(
  mu_values = seq(0, 0.8, by=0.2),
  k_values = c(5, 15, 30),
  tau_squared_values = c(0, 0.013333, 0.04, 0.12, 0.36),
  gamma = c(1.5),
  stringsAsFactors = FALSE
)

# Number of datasets simulated per parameter combination
n_sim <- 3200

# Initialize lists to store results
avg_trials_rep_list <- list()
all_datasets <- list()
param_values <- list()

# Set up parallel processing
num_cores <- 48
registerDoParallel(cores = num_cores)

# Loop over each parameter combination
for (i in 1:nrow(parDATA_MAR)) {
  
  cat("The", i, "th parameter combination out of", nrow(parDATA_MAR), "initiated.\n")
  
  tryCatch({
    
    # Extract parameter values
    mu <- parDATA_MAR$mu_values[i]
    k <- parDATA_MAR$k_values[i]
    tau_squared <- parDATA_MAR$tau_squared_values[i]
    gamma <- parDATA_MAR$gamma[i]
    
    # Initialize lists to store results for each parameter combination
    trials_list <- c()
    rep_list <- c()
    datasets <- list()
    
    # Generate datasets for each parameter combination
    for (j in 1:n_sim) {
      tryCatch({
        result <- simulate_ORB_data_MAR(n_studies = k,
                                    mu = mu,
                                    tau_squared = tau_squared,
                                    gamma = gamma)
        trials_list <- c(trials_list, result$Ntrials)
        rep_list <- c(rep_list, result$Nrep)
        
        # Add parameter values to dataset
        result$dataORB$mu_values <- mu
        result$dataORB$k_values <- k
        result$dataORB$tau_squared_values <- tau_squared
        result$dataORB$gamma <- gamma
        result$dataORB$average_Ntrials <- result$Ntrials
        result$dataORB$average_Nreported <- result$Nrep
        
        datasets[[j]] <- result$dataORB
      }, error = function(e) {
        cat("Error occurred at parameter combination", i, ": ", conditionMessage(e), "\n")
        stop("Error occurred, stopping execution")
      })
    }
    
    # Calculate average number of trials and repetitions for this parameter combination
    avg_trials <- mean(trials_list)
    avg_rep <- mean(rep_list)
    
    # Store averages for this parameter combination
    avg_trials_rep_list[[i]] <- c(mu, k, tau_squared, gamma, avg_trials, avg_rep)
    
    # Combine datasets for this parameter combination
    combined_datasets <- do.call(rbind, datasets)
    
    # Store combined datasets
    all_datasets[[i]] <- combined_datasets
    
    # Parallel computation for the current parameter combination
    res <- foreach(meta_data_miss = datasets, .options.RNG = seed, .combine = 'c') %dorng% {
      
      resultLL <- reORBadj(y = meta_data_miss$y_ORB,
                           s = meta_data_miss$s_ORB,
                           n1 = meta_data_miss$n_treatment,
                           n2 = meta_data_miss$n_control,
                           outcome = "benefit",
                           init_param = c(0.7, 0.37),
                           alpha_ben = 0.05,
                           alpha_ben_one.sided = TRUE, 
                           true.SE = NULL,
                           LR.CI = TRUE,
                           rho1 = 3,
                           rho2 = 3,
                           selection.benefit = "Constant.Constant",
                           opt_method = "L-BFGS-B")
      
      resultLC <- reORBadj(y = meta_data_miss$y_ORB,
                           s = meta_data_miss$s_ORB,
                           n1 = meta_data_miss$n_treatment,
                           n2 = meta_data_miss$n_control,
                           outcome = "benefit",
                           init_param = c(0.7, 0.37),
                           alpha_ben = 0.05,
                           alpha_ben_one.sided = TRUE,
                           true.SE = NULL,
                           LR.CI = TRUE,
                           rho1 = 3,
                           rho2 = 3,
                           selection.benefit = "Constant.Continous",
                           opt_method = "L-BFGS-B")
      
      resultCL <- reORBadj(y = meta_data_miss$y_ORB,
                           s = meta_data_miss$s_ORB,
                           n1 = meta_data_miss$n_treatment,
                           n2 = meta_data_miss$n_control,
                           outcome = "benefit",
                           init_param = c(0.7, 0.37),
                           alpha_ben = 0.05,
                           alpha_ben_one.sided = TRUE, 
                           true.SE = NULL,
                           LR.CI = TRUE,
                           rho1 = 3,
                           rho2 = 3,
                           selection.benefit = "Continous.Constant",
                           opt_method = "L-BFGS-B")
      
      resultCC1 <- reORBadj(y = meta_data_miss$y_ORB,
                            s = meta_data_miss$s_ORB,
                            n1 = meta_data_miss$n_treatment,
                            n2 = meta_data_miss$n_control,
                            outcome = "benefit",
                            init_param = c(0.7, 0.37),
                            alpha_ben = 0.05,
                            alpha_ben_one.sided = TRUE,
                            true.SE = NULL,
                            LR.CI = TRUE,
                            rho1 = 7,
                            rho2 = 1.5,
                            selection.benefit = "Continous.Continous",
                            opt_method = "L-BFGS-B")
      
      resultCC2 <- reORBadj(y = meta_data_miss$y_ORB,
                            s = meta_data_miss$s_ORB,
                            n1 = meta_data_miss$n_treatment,
                            n2 = meta_data_miss$n_control,
                            outcome = "benefit",
                            init_param = c(0.7, 0.37),
                            alpha_ben = 0.05,
                            alpha_ben_one.sided = TRUE,
                            true.SE = NULL,
                            LR.CI = TRUE,
                            rho1 = 1.5,
                            rho2 = 7,
                            selection.benefit = "Continous.Continous",
                            opt_method = "L-BFGS-B")
      
      
      
      list(LL = resultLL,
           LC = resultLC,
           CL = resultCL,
           CC1 = resultCC1,
           CC2 = resultCC2)
           
    }
    
    # Helper function to reshape the data from the parallelized loop
    formatting <- function(r, method) {
      new.list <- as.list(unlist(unname(r[names(r) == method])))
    }
    
    # Obtain lists for each adjustment method
    resultLL_list <- formatting(res, "LL")
    resultLC_list <- formatting(res, "LC")
    resultCL_list <- formatting(res, "CL")
    resultCC1_list <- formatting(res, "CC1")
    resultCC2_list <- formatting(res, "CC2")
    
    
    # Extract helper function
    extract <- function(lst, var.name) {
      l <- as.numeric(as.list(unlist(lst)[names(unlist(lst)) == var.name]))
      return(l)
    }
    
    # Extract theta_hat and CI for each method from the result lists 
    theta_hat_LL <- extract(resultLL_list, "mu_adjusted_benefit")
    theta_hat_low_LL <- extract(resultLL_list, "LR_mu_adjusted_low")
    theta_hat_up_LL <- extract(resultLL_list, "LR_mu_adjusted_up")
    tau2_hat_LL <- extract(resultLL_list, "tau_squared_adjusted")
    tau2_hat_low_LL <- extract(resultLL_list, "LR_tau_squared_adjusted_low")
    tau2_hat_up_LL <- extract(resultLL_list, "LR_tau_squared_adjusted_up")
    
    theta_hat_LC <- extract(resultLC_list, "mu_adjusted_benefit")
    theta_hat_low_LC <- extract(resultLC_list, "LR_mu_adjusted_low")
    theta_hat_up_LC <- extract(resultLC_list, "LR_mu_adjusted_up")
    tau2_hat_LC <- extract(resultLC_list, "tau_squared_adjusted")
    tau2_hat_low_LC <- extract(resultLC_list, "LR_tau_squared_adjusted_low")
    tau2_hat_up_LC <- extract(resultLC_list, "LR_tau_squared_adjusted_up")
    
    theta_hat_CL <- extract(resultCL_list, "mu_adjusted_benefit")
    theta_hat_low_CL <- extract(resultCL_list,  "LR_mu_adjusted_low")
    theta_hat_up_CL <- extract(resultCL_list,  "LR_mu_adjusted_up")
    tau2_hat_CL <- extract(resultCL_list, "tau_squared_adjusted")
    tau2_hat_low_CL <-extract(resultCL_list,  "LR_tau_squared_adjusted_low")
    tau2_hat_up_CL <- extract(resultCL_list, "LR_tau_squared_adjusted_up")
    
    theta_hat_CC1 <- extract(resultCC1_list, "mu_adjusted_benefit")
    theta_hat_low_CC1 <- extract(resultCC1_list,  "LR_mu_adjusted_low")
    theta_hat_up_CC1 <- extract(resultCC1_list,  "LR_mu_adjusted_up")
    tau2_hat_CC1 <- extract(resultCC1_list, "tau_squared_adjusted")
    tau2_hat_low_CC1 <-extract(resultCC1_list,  "LR_tau_squared_adjusted_low")
    tau2_hat_up_CC1 <- extract(resultCC1_list, "LR_tau_squared_adjusted_up")
    
    theta_hat_CC2 <- extract(resultCC2_list, "mu_adjusted_benefit")
    theta_hat_low_CC2 <- extract(resultCC2_list,  "LR_mu_adjusted_low")
    theta_hat_up_CC2 <- extract(resultCC2_list,  "LR_mu_adjusted_up")
    tau2_hat_CC2 <- extract(resultCC2_list, "tau_squared_adjusted")
    tau2_hat_low_CC2 <-extract(resultCC2_list,  "LR_tau_squared_adjusted_low")
    tau2_hat_up_CC2 <- extract(resultCC2_list, "LR_tau_squared_adjusted_up")
    
    
    theta_hat_Unadj <- extract(resultLL_list,  "mu_unadjusted")
    theta_hat_low_Unadj <- extract(resultLL_list, "LR_mu_unadjusted_low")
    theta_hat_up_Unadj <- extract(resultLL_list, "LR_mu_unadjusted_up")
    tau2_hat_Unadj <- extract(resultLL_list, "tau_squared_adjusted")
    tau2_hat_low_Unadj <-extract(resultLL_list,  "LR_tau_squared_adjusted_low")
    tau2_hat_up_Unadj <- extract(resultLL_list, "LR_tau_squared_adjusted_up")
    
    
    # Calculate the summary statistics for each method
    summaryLL <- calculate_summary(theta_hat_LL, theta_hat_low_LL, theta_hat_up_LL, mu)
    summaryLC <- calculate_summary(theta_hat_LC, theta_hat_low_LC, theta_hat_up_LC, mu)
    summaryCL <- calculate_summary(theta_hat_CL, theta_hat_low_CL, theta_hat_up_CL, mu)
    summaryCC1 <- calculate_summary(theta_hat_CC1, theta_hat_low_CC1, theta_hat_up_CC1, mu)
    summaryCC2 <- calculate_summary(theta_hat_CC2, theta_hat_low_CC2, theta_hat_up_CC2, mu)
    summaryUnadj <- calculate_summary(theta_hat_Unadj, theta_hat_low_Unadj, theta_hat_up_Unadj, mu)
    
    # Store average results in parDATA_MAR dataframe
    parDATA_MAR[i, paste0("LL_", names(summaryLL))] <- unlist(summaryLL)
    parDATA_MAR[i, paste0("LC_", names(summaryLC))] <- unlist(summaryLC)
    parDATA_MAR[i, paste0("CL_", names(summaryCL))] <- unlist(summaryCL)
    parDATA_MAR[i, paste0("CC1_", names(summaryCC1))] <- unlist(summaryCC1)
    parDATA_MAR[i, paste0("CC2_", names(summaryCC2))] <- unlist(summaryCC2)
    parDATA_MAR[i, paste0("Unadj_", names(summaryUnadj))] <- unlist(summaryUnadj)
    
    summaryLL_tau2 <- calculate_summary(tau2_hat_LL, tau2_hat_low_LL, tau2_hat_up_LL, tau_squared)
    summaryLC_tau2 <- calculate_summary(tau2_hat_LC, tau2_hat_low_LC, tau2_hat_up_LC, tau_squared)
    summaryCL_tau2 <- calculate_summary(tau2_hat_CL, tau2_hat_low_CL, tau2_hat_up_CL, tau_squared)
    summaryCC1_tau2 <- calculate_summary(tau2_hat_CC1, tau2_hat_low_CC1, tau2_hat_up_CC1, tau_squared)
    summaryCC2_tau2 <- calculate_summary(tau2_hat_CC2, tau2_hat_low_CC2, tau2_hat_up_CC2, tau_squared)
    summaryUnadj_tau2 <- calculate_summary(tau2_hat_Unadj, tau2_hat_low_Unadj, tau2_hat_up_Unadj, tau_squared)
    
    # Store average results in parDATA_MAR dataframe
    parDATA_MAR[i, paste0("tau2_LL_", names(summaryLL_tau2))] <- unlist(summaryLL_tau2)
    parDATA_MAR[i, paste0("tau2_LC_", names(summaryLC_tau2))] <- unlist(summaryLC_tau2)
    parDATA_MAR[i, paste0("tau2_CL_", names(summaryCL_tau2))] <- unlist(summaryCL_tau2)
    parDATA_MAR[i, paste0("tau2_CC1_", names(summaryCC1_tau2))] <- unlist(summaryCC1_tau2)
    parDATA_MAR[i, paste0("tau2_CC2_", names(summaryCC2_tau2))] <- unlist(summaryCC2_tau2)
    parDATA_MAR[i, paste0("tau2_Unadj_", names(summaryUnadj_tau2))] <- unlist(summaryUnadj_tau2)
    
    
    # Clear intermediate objects and free memory
    rm(
      resultLL_list, 
      resultLC_list, 
      resultCL_list, 
      resultCC1_list, 
      resultCC2_list)
     
    gc()
    
    
    # Save parDATA_MAR until this point
    saveRDS(parDATA_MAR, file = paste0("parDATA_MAR15_", i, ".rds"))
    
    
    # Remove previous intermediate files
    if (i > 1) {
      previous_files <- paste0("parDATA_MAR15_", 1:(i-1), ".rds")
      file.remove(previous_files)
    }
    
    # Clear intermediate datasets and free memory
    rm(datasets)
    gc()
    
    
  }, error = function(e) {
    cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
    break
  })
  
}

# If any error occurred, return parDATA_MAR until that point
if (exists("i") && i < nrow(parDATA_MAR)) {
  parDATA_MAR <- parDATA_MAR[1:i, ]
}

# Save results
saveRDS(parDATA_MAR, file = "AdjustResults3200k_15_MAR.rds")


