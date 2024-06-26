
#The function reORBadj adjusts for ORB according to the selection model framework developed in this project
#It can be seen as an extension of the Copas et al. (2019) methodology, for the random effects model, and without
#the need for an ORBIT classification of unreported study outcomes.

#Inputs
#a - vector for observed counts in the treatment arm for each study in the meta-analysis, unreported study outcomes have entry "unrep"
#c - vector for observed counts in the control arm for each study in the meta-analysis, unreported study outcomes have entry "unrep"
#mu1 - vector for observed means in the treatment arm for each study in the meta-analysis, unreported study outcomes have entry "unrep"
#mu2 - vector for observed means in the control arm for each study in the meta-analysis, unreported study outcomes have entry "unrep"
#sd1 - vector for observed standard errors in the treatment arm for each study in the meta-analysis, unreported study outcomes have entry "unrep"
#sd2 - vector for observed standard errors in the control arm for each study in the meta-analysis, unreported study outcomes have entry "unrep"
#y - vector for observed treatment effects (normally distributed, eg log RR) for each study in the meta-analysis, unreported study outcomes have entry "unrep"
#s - vector for observed standard errors of the the treatment effect for each study in the meta-analysis, unreported study outcomes have entry "unrep"
#n1 - vector for sample sizes in treatment arm, must be reported for all studies
#n2 - vector for sample sizes in control arm, must be reported for all studies
#outcome - string character "beneficial"
#init_param - vector of length two, containing the values to be used for initialization of the joint optimization of the treatment effect and heterogeneity variance
#alpha_ben - value for desired confidence level, e.g., 0.05
#alpha_ben_one.sided - if TRUE, one sided significance is used, FALSE for two-sided
#true.SE - if the true standard errors are known, valid only for the case in which we have y and s as inputs; in this case the SE are not imputed
#LR.CI - if TRUE, calculates the profile likelihood confidence intervals using the likelihood ratio statistic
#Wald.CI - if TRUE, calculates Wald CI 
#selection.benefit - string factor specifying the selection function to be used for the ORB adjustment 
#                  "Constant-constant" W_A in manuscript
#                  "Constant-continous" W_B in manuscript
#                  "Continous-constant" W_C in manuscript
#                 "Continous-continous" W_D in manuscript
#                  "DGM" W_DGM in manuscript
#rho1 - parameter for continuous part of selection function W_B in manuscript, "Constant-continous"
#rho2 - parameter for continuous part of selection function W_C in manuscript, "Continous-constant"
#opt_method - string character for the desired optimizatin method, e.g., "L-BFGS-B"
#lower - vector of length 2 for the lower bound of the optimization if "L-BFGS-B" is used
#upper - vector of length 2 for the upper bound of the optimization if "L-BFGS-B" is used
#gam - parameter used in DGM and consequently to be used in DGM selection function for correct model specification

reORBadj <- function(a=NULL, 
                     c=NULL,
                     mu1=NULL, 
                     mu2=NULL, 
                     sd1=NULL, 
                     sd2=NULL,
                     y=NULL, 
                     s=NULL,
                     n1,
                     n2,
                     outcome,
                     init_param,
                     alpha_ben=NULL,
                     alpha_ben_one.sided = TRUE,
                     true.SE=NULL, 
                     LR.CI = FALSE,
                     Wald.CI = FALSE,
                     selection.benefit = "Copas.oneside",
                     rho1=3,
                     rho2=3,
                     opt_method="L-BFGS-B",
                     lower = c(-5, 0.00001),
                     upper = c(5,5),
                     gam = 1.5
) {
  
  #Can take in binary counts and calculate log RR: a, c
  #Means and calculate differences in means: y1, y2, sd1, sd2
  #Directly normally distributed effect measure: y, s
  
  #Beneficial outcome, i.e., positive direction of treatment
  
  #Selection functions for the probability of reporting as a function of pvalue
  #SS: Constant-constant
  #ST: Constant-continous
  #TS: Continous-constant
  #TT: Continous-continous
  #DGM: Corresponds to the function to use to generate ORB
  
  #If alpha_ben_one.sided=TRUE, we use the one-sided threshold 
  sel.ben <- selection.benefit
  
  #Parameters for continuous parts of selection functions
  rho1 <- rho1
  rho2 <- rho2
  
  #Parameter used in DGM and consequently to be used in DGM selection function for correct model specification
  gam <- gam
  
  #If FALSE, we do not calculate the Wald CI - they are numerically unstable sometimes
  if (Wald.CI){
    my.hessian <- TRUE
  } else{
    my.hessian <- FALSE
  }
  
  #Optimization method
  method <- opt_method
  lower <- lower
  upper <- upper
  
  
  #Binary input data to calculate log RR
  if (!is.null(a) & !is.null(c)){
    
    #Reported outcomes
    #Indecies where we do not have unrepoted outcomes, i.e., we have reported outcomes
    #If we turn C into numeric, the unreported become NA
    
    Rep_index <- which(!is.na(as.numeric(a)))
    HR_index <- which(a == "unrep")
    
    # a,c,n1,n2, n values for the reported studies
    a_rep <- as.numeric(a[Rep_index])
    c_rep <- as.numeric(c[Rep_index])
    n1_rep <- as.numeric(n1[Rep_index])
    n2_rep <- as.numeric(n2[Rep_index])
    
    if (length(a_rep == 0) | length(c_rep == 0) > 0){ #Continuity correction
      
      a_0_index <- c(which(a_rep == 0), which(c_rep == 0)) #Where the zero cell counts are
      a_rep[a_0_index] <- a_rep[a_0_index] + 0.5
      c_rep[a_0_index] <- c_rep[a_0_index] + 0.5
      n1_rep[a_0_index] <- n1_rep[a_0_index] + 0.5
      n2_rep[a_0_index] <- n2_rep[a_0_index] + 0.5 #Add 0.5 to the cell counts with 0
      
      a_0_both <- which(a_rep == 0.5 & c_rep == 0.5)
      
      if (length(a_0_both)>0){
        
        a_rep <- a_rep[-a_0_both]
        c_rep <- c_rep[-a_0_both]
        n1_rep <- n1_rep[-a_0_both]
        n2_rep <- n2_rep[-a_0_both]
        
      } else {
        
        a_rep <- a_rep
        c_rep <- c_rep
        n1_rep <- n1_rep
        n2_rep <- n2_rep
        
      }
      
    } else {
      
      a_rep <- a_rep
      c_rep <- c_rep
      n1_rep <- n1_rep
      n2_rep <- n2_rep
    }
    
    
    #How many studies are reported?
    N_rep <- length(Rep_index)
    
    #Unreported study sizes, we might have info from n1,n2 or just the total
    ntot <- as.numeric(n1) + as.numeric(n2)
    n_HR <- as.numeric(ntot[HR_index])
    
    
    #log RR and standard error
    logRR <- log((a_rep*n2_rep)/(c_rep*n1_rep))
    s <- sqrt(((n1_rep - a_rep)/(n1_rep*a_rep)) + ((n2_rep - c_rep)/(n2_rep*c_rep)))
    sigma_squared <- s^2
    
    
    
    #Imputed values of sigma squared for the unreported studies
    #K value based on the reported studies, see Copas et al. (2019)
    k <- sum(1/sigma_squared)/sum((n1_rep + n2_rep))
    
    #Mean differences as input data (continuous)
    #Same procedure
    
  } else if (!is.null(mu1) & !is.null(mu2) & !is.null(sd1) & !is.null(sd2)){
    
    Rep_index <- which(!is.na(as.numeric(mu1)))
    
    HR_index <- which(mu1 == "unrep")
      
    #Unreported study sizes, we might have info from n1,n2 or just the total
    ntot <- as.numeric(n1) +as.numeric(n2)
    n_HR <- as.numeric(ntot[HR_index])
    
    #mu1,mu2,n1,n2,n values for the reported studies
    mu1_rep <- as.numeric(mu1[Rep_index])
    mu2_rep <- as.numeric(mu2[Rep_index])
    sd1_rep <- as.numeric(sd1[Rep_index])
    sd2_rep <- as.numeric(sd2[Rep_index])
    n1_rep <- as.numeric(n1[Rep_index])
    n2_rep <- as.numeric(n2[Rep_index])
    
    #How many studies are reported?
    N_rep <- length(Rep_index)
    
    #Differenece in means
    #Standard errors given
    logRR <- mu1_rep - mu2_rep
    s <- sqrt((as.numeric(sd1_rep)^2)/(as.numeric(n1_rep)) + (as.numeric(sd2_rep)^2)/(as.numeric(n2_rep)))
    sigma_squared <- s^2
    
    k <- sum(1/sigma_squared)/sum((n1_rep + n2_rep))
    
    #We directly have the treatment effect and standard error
    
  } else if (!is.null(y) & !is.null(s)) {
    
    #Indecies where we have the reported outcomes and the unreported HR
    Rep_index <- which(!is.na(as.numeric(y)))
    
    HR_index <- which(y == "unrep")
    
    #Observed treatment effect, standard error^2, and sample sizes of reported outcomes
    logRR <- as.numeric(y[Rep_index])
    sigma_squared <- (as.numeric(s[Rep_index]))^2
    n1_rep <- as.numeric(n1[Rep_index])
    n2_rep <- as.numeric(n2[Rep_index])
    
    
    #Total sample size of unreported outcomes
    n_HR <- as.numeric(n1[HR_index]) + as.numeric(n2[HR_index])
    
    
    #k estimation, based on reported studies
    k <- sum(1/sigma_squared)/sum((n1_rep + n2_rep))
    
    
    
  } else {
    
    return("Error: invalid inputs. Input either a and c values or y1,sd1 and y2,sd2 values.")
  }
  
  
  #Possibility to pass the true SE to the function
  if (!is.null(true.SE)){
    
    sigma_squared_imputed <- (as.numeric(true.SE)[HR_index])^2
  } else {
    #Imputed variances for the HR studies
    #See Copas et al. (2019)
    sigma_squared_imputed <- 1/(k*n_HR)
    
  }
  
  #Average sigma squared value that is used when we do not adjust for ORB and when we do
  sigma_squared_average_unadjusted <- mean(sigma_squared)
  sigma_squared_average_adjusted <- mean(c(sigma_squared, sigma_squared_imputed))
  
  #Unadjusted log-likelihood function to be maximized
  f.u <- function(params, logRR, sigma_squared) {
    mu <- params[1]
    tau_squared <- params[2]
    if (tau_squared < 0 ){
      - Inf
    } else {
      
      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
    }
  }
  
  #Set initial values for mu and tau_squared
  init_params <- init_param
  
  #Maximize unadjusted function optim()
  if (method == "L-BFGS-B"){
    
    fit.u <- optim(init_params, f.u, logRR = logRR, sigma_squared = sigma_squared,
                   #method = "Nelder-Mead",
                   method=method,
                   lower = lower,
                   upper = upper,
                   control = list(fnscale = -1),
                   hessian=my.hessian)
  } else {
    
    fit.u <- optim(init_params, f.u, logRR = logRR, sigma_squared = sigma_squared,
                   method = method,
                   control = list(fnscale = -1),
                   hessian=my.hessian)
    
  }
  
  #Return unadjusted mu and tau_squared
  mle.u <- fit.u$par[1]
  mle.tau <- max(fit.u$par[2],0)
  
  
  #Beneficial outcome adjustment for ORB
  
  if(outcome == "benefit"){
    
    z_alpha.copas <- qnorm(1-alpha_ben/2) #two-sided threshold
    
    if (alpha_ben_one.sided == TRUE){ #one-sided threshold (reccomended for consistency with simulation process and literature)
      
      z_alpha <- qnorm(1-alpha_ben)
    } else {
      z_alpha <- qnorm(1-alpha_ben/2)
    }
    
    #Adjusted log-likelihood function for beneficial outcome to be maximized
    f.adj.b <- function(params, logRR, sigma_squared, sigma_squared_imputed) {
      mu <- params[1]
      tau_squared <- params[2]
      
      if(tau_squared < 0){
        - Inf
      } else {
        
        #The contribution from the reported studies is always present
        -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +
          
          
          if (length(sigma_squared_imputed) > 0) {
            
            if (sel.ben == "Constant.Constant"){
              
              #sum(log(pnorm((z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))))
              
              #Probability of reporting
              SS <- function(y, sigma_squared_imputed, alpha_ben) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                ifelse(p <= alpha_ben, 
                       1,
                       0)
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, alpha_ben) {
                
                weight <- 1 - SS(y, sigma_squared_imputed, alpha_ben)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, alpha_ben), #change here April 2. added Vectorize
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
              
            } else if (sel.ben == "Constant.Continous") {
              
              #Probability of reporting
              TS <- function(y, sigma_squared_imputed, alpha_ben) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                ifelse(p <= alpha_ben, 
                       1,
                       (p^(-rho2))/(alpha_ben^(-rho2)))
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, alpha_ben) {
                
                weight <- 1 - TS(y, sigma_squared_imputed, alpha_ben)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, alpha_ben), #change here April 2. added Vectorize
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
            } else if (sel.ben == "Continous.Constant") {
              
              #Probability of reporting
              ST <- function(y, sigma_squared_imputed, alpha_ben) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                ifelse(p <= alpha_ben, 
                       1-(p^rho1)/(alpha_ben^rho1),
                       0
                )
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, alpha_ben) {
                
                weight <- 1 - ST(y, sigma_squared_imputed, alpha_ben)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, alpha_ben),
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
            } else if (sel.ben == "Continous.Continous") {
              
              #Probability of reporting
              TT <- function(y, sigma_squared_imputed, alpha_ben) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                ifelse(p <= alpha_ben, 
                       1 - (1-0.5)*(p^rho1)/(alpha_ben^rho1),
                       0.5*(p^(-rho2))/(alpha_ben^(-rho2))
                )
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, alpha_ben) {
                
                weight <- 1 - TT(y, sigma_squared_imputed, alpha_ben)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, alpha_ben),
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
            
            } else if (sel.ben == "DGM") {
              
              #Probability of reporting
              DGM <- function(y, sigma_squared_imputed, gam) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                exp(-4*p^(gam))
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, gam) {
                
                weight <- 1 - DGM(y, sigma_squared_imputed, gam)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, gam),
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
              
            } else {
              
              stop("Error: Invalid weight function specified. ")
            }
            
            
            
          } else {
            0  # Return 0 if sigma_squared_imputed is empty
          }
        
      }
      
      
    }
    
    #Maximize log-likelihood
    if (method == "L-BFGS-B"){
      
      fit.adj.b <- optim(init_params, f.adj.b, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                         method = method,
                         lower = lower,
                         upper = upper,
                         control = list(fnscale = -1),
                         
                         hessian=my.hessian)
    } else {
      
      fit.adj.b <- optim(init_params, f.adj.b, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                         method = method,
                         control = list(fnscale = -1),
                         hessian=my.hessian)
      
      
    }
    
    #Return adjusted mu and tau_squared
    mle.b <- fit.adj.b$par[1]
    mle.b.tau <- max(fit.adj.b$par[2],0)
    
    
    #LIKELIHOOD RATIO CONFIDENCE INTERVALS
    
    if (LR.CI){
      
      #Unadjusted
      z <- qchisq(1-alpha_ben, df=1) #3.841
      
      #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
      ll.u <- function(mu, tau_squared, logRR, sigma_squared) {
        
        -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
      }
      
      #Profile log-likelihood
      pl.u <- function(mu, logRR, sigma_squared) { #take in vector of mus
        
        res <- mu
        
        for (i in seq_along(mu)) { #for all these values of mu
          optimResult <- optim(par = init_param[2],
                               fn = function(tau_squared) ll.u(mu[i], tau_squared, logRR=logRR, sigma_squared = sigma_squared),
                               method = "Brent",
                               lower= 0.0001,
                               upper=10,
                               control = list(fnscale = -1))
          
          res[i] <- optimResult$value
        }
        return(res)
      }
      
      f <- function(mu, logRR, sigma_squared){
        pl.u(mu, logRR=logRR, sigma_squared=sigma_squared) - pl.u(mle.u, logRR=logRR, sigma_squared=sigma_squared) + 1/2*qchisq(0.95, df=1)
      }
      
      
      eps <- sqrt(.Machine$double.eps)
      lowerBound.u <- uniroot(f, interval = c(-5, mle.u), logRR=logRR, sigma_squared=sigma_squared)$root
      upperBound.u <- uniroot(f, interval = c( mle.u, 5), logRR=logRR, sigma_squared=sigma_squared)$root
      
      #Profile log-likelihood for tau squared
      pl.u.tau <- function(tau_squared, logRR, sigma_squared) { #take in vector of tau_squares
        
        res <- tau_squared
        
        for (i in seq_along(tau_squared)) { #for all these values of mu
          optimResult <- optim(par = init_param[1],
                               fn = function(mu) ll.u(tau_squared[i], mu, logRR=logRR, sigma_squared = sigma_squared),
                               method = "Brent",
                               lower= -5,
                               upper=5,
                               control = list(fnscale = -1))
          
          res[i] <- optimResult$value
        }
        return(res)
      }
      
      f.tau <- function(tau_squared, logRR, sigma_squared){
        pl.u.tau(tau_squared, logRR=logRR, sigma_squared=sigma_squared) - pl.u.tau(mle.tau, logRR=logRR, sigma_squared=sigma_squared) + 1/2*qchisq(0.95, df=1)
      }
      
      
      eps <- sqrt(.Machine$double.eps)
      lowerBound.u.tau <- max(uniroot(f.tau, interval = c(-10, mle.tau), logRR=logRR, sigma_squared=sigma_squared)$root,0)
      upperBound.u.tau <- max(uniroot(f.tau, interval = c( mle.tau, 10), logRR=logRR, sigma_squared=sigma_squared)$root,0)
      
      
      #Adjusted benefit
      
      #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
      ll.b <- function(mu, tau_squared, logRR, sigma_squared, sigma_squared_imputed) {
        
        if(tau_squared < 0){
          - Inf
        } else {
          
        #The contribution from the reported studies is always present
        -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +
          
          
          if (length(sigma_squared_imputed) > 0) {
            
            if (sel.ben == "Constant.Constant"){
              
              #sum(log(pnorm((z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))))
              
              #Probability of reporting
              SS <- function(y, sigma_squared_imputed, alpha_ben) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                ifelse(p <= alpha_ben, 
                       1,
                       0)
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, alpha_ben) {
                
                weight <- 1 - SS(y, sigma_squared_imputed, alpha_ben)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, alpha_ben), #change here April 2. added Vectorize
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
              
            } else if (sel.ben == "Constant.Continous") {
              
              #Probability of reporting
              TS <- function(y, sigma_squared_imputed, alpha_ben) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                ifelse(p <= alpha_ben, 
                       1,
                       (p^(-rho2))/(alpha_ben^(-rho2)))
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, alpha_ben) {
                
                weight <- 1 - TS(y, sigma_squared_imputed, alpha_ben)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, alpha_ben), #change here April 2. added Vectorize
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
            } else if (sel.ben == "Continous.Constant") {
              
              #Probability of reporting
              ST <- function(y, sigma_squared_imputed, alpha_ben) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                ifelse(p <= alpha_ben, 
                       1-(p^rho1)/(alpha_ben^rho1),
                       0
                )
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, alpha_ben) {
                
                weight <- 1 - ST(y, sigma_squared_imputed, alpha_ben)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, alpha_ben),
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
            } else if (sel.ben == "Continous.Continous") {
              
              #Probability of reporting
              TT <- function(y, sigma_squared_imputed, alpha_ben) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                ifelse(p <= alpha_ben, 
                       1 - (1-0.5)*(p^rho1)/(alpha_ben^rho1),
                       0.5*(p^(-rho2))/(alpha_ben^(-rho2))
                )
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, alpha_ben) {
                
                weight <- 1 - TT(y, sigma_squared_imputed, alpha_ben)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, alpha_ben),
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
              
              
              
            } else if (sel.ben == "DGM") {
              
              #Probability of reporting
              DGM <- function(y, sigma_squared_imputed, gam) {
                
                p = pnorm(-y/sqrt(sigma_squared_imputed))
                
                exp(-4*p^(gam))
              }
              
              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, gam) {
                
                weight <- 1 - DGM(y, sigma_squared_imputed, gam)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }
              
              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, gam),
                          lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
              })))
              
              
            } else {
              
              stop("Error: Invalid weight function specified. ")
            }
            
            
            
          } else {
            0  # Return 0 if sigma_squared_imputed is empty
          }
          
        }
        
        
      }
      
      #Adjusted profile log likelihood for mu
      pl.b <- function(mu, logRR, sigma_squared, sigma_squared_imputed) { #take in vector of mus
        
        res <- mu
        
        for (i in seq_along(mu)) { #for all these values of mu
          optimResult <- optim(par = init_param[2],
                               fn = function(tau_squared) ll.b(mu[i], tau_squared,
                                                               logRR=logRR,
                                                               sigma_squared = sigma_squared,
                                                               sigma_squared_imputed = sigma_squared_imputed),
                               method = "Brent",
                               lower=0.00001,
                               upper=10,
                               
                               control = list(fnscale = -1))
          
          res[i] <- optimResult$value
        }
        return(res)
      }
      
      f.b <- function(mu, logRR, sigma_squared, sigma_squared_imputed){
        
        pl.b(mu, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) - pl.b(mle.b, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) + 1/2*qchisq(0.95, df=1)
      }
      
      lowerBound.b <- uniroot(f.b, interval = c(-10, mle.b), logRR=logRR,
                              sigma_squared=sigma_squared,
                              sigma_squared_imputed=sigma_squared_imputed)$root
      upperBound.b <- uniroot(f.b, interval = c(mle.b, 10), logRR=logRR,
                              sigma_squared=sigma_squared,
                              sigma_squared_imputed = sigma_squared_imputed)$root
      
      #Adjusted profile log likelihood for tau squared
      pl.b.tau <- function(tau_squared, logRR, sigma_squared, sigma_squared_imputed) { #take in vector of tau_squares
        
        res <- tau_squared
        
        for (i in seq_along(tau_squared)) { #for all these values of tau_squared
          optimResult <- optim(par = init_param[1],
                               fn = function(mu) ll.b(tau_squared[i], mu,
                                                               logRR=logRR,
                                                               sigma_squared = sigma_squared,
                                                               sigma_squared_imputed = sigma_squared_imputed),
                               method = "Brent",
                               lower=-5,
                               upper=5,
                               
                               control = list(fnscale = -1))
          
          res[i] <- optimResult$value
        }
        return(res)
      }
      
      f.b.tau <- function(tau_squared, logRR, sigma_squared, sigma_squared_imputed){
        
        pl.b.tau(tau_squared, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) - pl.b.tau(mle.b.tau, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) + 1/2*qchisq(0.95, df=1)
      }
      
      lowerBound.b.tau <- 0
      
      tryCatch({
        
      lowerBound.b.tau <- 
        uniroot(f.b.tau, interval = c(-5, fit.adj.b$par[2]), logRR=logRR,
                              sigma_squared=sigma_squared, extendInt = "yes")$root
      }, error = function(e) {
        
        lowerBound.b.tau <- 0})
             
      upperBound.b.tau <- 0
      
      tryCatch({                        
      
      upperBound.b.tau <- 
        uniroot(f.b.tau, interval = c(fit.adj.b$par[2], 5), logRR=logRR,
                              sigma_squared=sigma_squared, sigma_squared_imputed=sigma_squared_imputed, extendInt="yes")$root
      }, error = function(e) {
        
        upperBound.b.tau <- 0})
   
      
      
      
      if (Wald.CI){
        
        
        #WALD CONFIDENCE INTERVALS
        a <- alpha_ben #for harm Copas et al use 99% conf level
        #Unadjusted
        fisher_info.u <- solve(-fit.u$hessian)
        s.u <- sqrt(diag(fisher_info.u)[1])
        ci.u <- fit.u$par[1] + qnorm(c(a/2, 1-a/2)) * s.u
        #Adjusted benefit
        fisher_info.adj.b <- solve(-fit.adj.b$hessian)
        s.adj.b <- sqrt(diag(fisher_info.adj.b)[1])
        ci.u.adj.b <- fit.adj.b$par[1] + qnorm(c(a/2, 1-a/2)) * s.adj.b
        
        
        
        return(list(mu_unadjusted = mle.u,
                    LR_mu_unadjusted_low = lowerBound.u,
                    LR_mu_unadjusted_up = upperBound.u,
                    
                    
                    CI_unadjusted_low_WALD = ci.u[1],
                    CI_unadjusted_up_WALD = ci.u[2],
                    
                    mu_adjusted_benefit = mle.b,
                    LR_mu_adjusted_low = lowerBound.b,
                    LR_mu_adjusted_up = upperBound.b,
                    
                    tau_squared_unadjusted = mle.tau,
                    LR_tau_squared_unadjusted_low = lowerBound.u.tau,
                    LR_tau_squared_unadjusted_up = upperBound.u.tau,
                    
                    tau_squared_adjusted = mle.b.tau,
                    LR_tau_squared_adjusted_low = lowerBound.b.tau,
                    LR_tau_squared_adjusted_up = upperBound.b.tau,
                    
                    
                    
                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted,
                    
                    
                    CI_adjusted_benefit_low_WALD = ci.u.adj.b[1],
                    CI_adjusted_benefit_up_WALD = ci.u.adj.b[2]
                    
        ))
        
      } else {
        
        return(list(mu_unadjusted = mle.u,
                    LR_mu_unadjusted_low = lowerBound.u,
                    LR_mu_unadjusted_up = upperBound.u,
                    
                    
                    
                    mu_adjusted_benefit = mle.b,
                    LR_mu_adjusted_low = lowerBound.b,
                    LR_mu_adjusted_up = upperBound.b,
                    
                    tau_squared_unadjusted = mle.tau,
                    LR_tau_squared_unadjusted_low = lowerBound.u.tau,
                    LR_tau_squared_unadjusted_up = upperBound.u.tau,
                    
                    tau_squared_adjusted = mle.b.tau,
                    LR_tau_squared_adjusted_low = lowerBound.b.tau,
                    LR_tau_squared_adjusted_up = upperBound.b.tau,
                    
                    
                    
                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted
                    
                    
                    
        ))
        
        
      }
      
      
    } else {
      
      
      
      return(list(
        mu_unadjusted = mle.u,
        mu_adjusted_benefit = mle.b,
        tau_squared_unadjusted = mle.tau,
        tau_squared_adjusted = mle.b.tau,
        average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
        sigma_squared_average_adjusted = sigma_squared_average_adjusted
      ))
      
      
    }
    
    
  } else {
    
    return("invalid outcome input")
  }
  
  
}

