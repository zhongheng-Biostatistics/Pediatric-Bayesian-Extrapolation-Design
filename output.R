#' output function
#' 
#' This function computes Bayesian posterior decision values based on MCMC samples 
#' obtained from `posterior_sample`. It supports both "type I" error and "power" analysis.
#' 
#' @param eps_H Threshold in H_0.
#' @param lower Lower bound of the region of interest.
#' @param upper Upper bound of the region of interest.
#' @param w_1 Vector representing weights for adult information in the REPP framework.
#' @param N Number of iterations for sampling.
#' @param beta_adu Adult's coefficients.
#' @param post_eps Vector of posterior probability thresholds.
#' @param posterior_sample MCMC sample matrix, output from `posterior_sample` function.
#' @param n.burnin Number of burn-in samples for JAGS.
#' @param n.iter Total number of MCMC iterations for JAGS.
#' @param n.thin Thinning factor for JAGS.
#' @param type Indicator specifying whether to compute "type I" error or "power".
#' 
#' @return If `type = "type I"`, returns a list containing the `post` matrix and an identifier.
#'         If `type = "power"`, returns the `post` matrix directly.
#' 

output <- function(eps_H,lower, upper, w_1, N, beta_adu, post_eps, posterior_sample, 
                   n.burnin = 10000, n.iter = 80000, n.thin = 20, type) {
  
  # Compute the number of rows per MCMC iteration
  n.row <- (n.iter - n.burnin) / n.thin
  
  # Initialize the Bayes decision matrix
  post <- matrix(0, ncol = length(post_eps), nrow = length(w_1))
  
  # Assign column names based on post_eps values
  colnames(post) <- paste0("max_post_", post_eps)
  
  # Assign row names based on w_1 values
  rownames(post) <- paste0("w=", w_1)
  
  # Define the function for computing maximum log odds difference
  comp_max <- function(beta1, beta2, lower, upper) {
    f <- function(x) {
      exp(beta1[1] + x * beta1[2]) / (1 + exp(beta1[1] + x * beta1[2]))
    }
    g <- function(x) {
      exp(beta2[1] + x * beta2[2]) / (1 + exp(beta2[1] + x * beta2[2]))
    }
    h <- function(x) { f(x) - g(x) }
    return(max(c((optimize(h, lower = lower, upper = upper, maximum = TRUE)$objective), h(lower), h(upper))))
  }
  
  # Compute posterior decision values
  for (l in 1:length(post_eps)) {
    print(paste0("l=", l))
    
    Bayes_post <- matrix(0, ncol = length(w_1), nrow = N)
    
    for (j in 1:length(w_1)) {
      for (i in 1:N) {
        count_max <- 0
        
        # Extract MCMC samples for the current (j, i) combination
        result <- posterior_sample[((j - 1) * n.row * N + n.row * (i - 1) + 1):((j - 1) * n.row * N + n.row * i), ]
        Pmatrix <- result[, 1:2]
        
        for (k in 1:nrow(Pmatrix)) {
          if (comp_max(beta_adu, Pmatrix[k, ], lower, upper) <= eps_H) {
            count_max <- count_max + 1
          }
        }
        
        Bayes_post[i, j] <- as.numeric(count_max / nrow(Pmatrix) > post_eps[l])
      }
    }
    
    post[, l] <- apply(Bayes_post, 2, mean)
  }
  
  # Return output based on the type
  if (type == "type I") {
    return(list(typeI = post))
  } else {
    return(list(power = post))
  }
}


output_results_typeI <- output(eps_H=0.2,
  lower = 2.5, upper = 5, w_1 = seq(0.1, 0.5, length = 5), 
  N =1, beta_adu = c(-2.830380, 1.412872),
   post_eps = c(0.8, 0.85, 0.9, 0.95, 0.99), 
  posterior_sample = MCMC_results_typeI,
   n.burnin = 10000, n.iter = 80000, n.thin = 20,
   type = "type I"
 )

output_results_power <- output(eps_H=0.2,
                         lower = 2.5, upper = 5, w_1 = seq(0.1, 0.5, length = 5), 
                         N =1, beta_adu = c(-2.830380, 1.412872),
                         post_eps = c(0.8, 0.85, 0.9, 0.95, 0.99), 
                         posterior_sample = MCMC_results_power,
                         n.burnin = 10000, n.iter = 80000, n.thin = 20,
                         type = "power"
)
