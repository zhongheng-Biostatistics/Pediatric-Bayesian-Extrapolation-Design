#' posterior_sample function
#' 
#' This function performs MCMC sampling using JAGS to estimate pediatric model coefficients.
#' It supports two modes: "type I" error and "power" analysis.
#' 
#' @param eps_H Threshold in H_0.
#' @param lower Lower bound of the region of interest.
#' @param upper Upper bound of the region of interest.
#' @param n_ped Pediatric sample size.
#' @param beta_adu Adult's coefficients.
#' @param type Indicator for the type of result, must be either "type I" or "power".
#' @param weight_beta Weights for different beta vectors.
#' @param weight_delta Only applicable when type = "power". Specifies weights for delta = eps_H * k / 10, k = 0,...,9.
#' @param beta_ped_all A collection of pediatric beta coefficients.
#' @param w_1 Weights for adult information in the REPP framework.
#' @param x_lower Lower bound of the exposure range.
#' @param x_upper Upper bound of the exposure range.
#' @param beta1_w Mixture normal Weight for beta1.
#' @param beta2_w Mixture normal Weight for beta2.
#' @param mu1_1 The first element of Mean vector for Mixture normal approximation of beta1.
#' @param prec1_1 The first element of precision vector for Mixture normal approximation of beta1.
#' @param mu1_2 The second element of mean for Mixture normal approximation of beta1.
#' @param prec1_2 The second element of precision vector for Mixture normal approximation of beta1.
#' @param mu2_1 The first element of Mean vector for Mixture normal approximation of beta2.
#' @param prec2_1 The first element of precision vector for Mixture normal approximation of beta2
#' @param mu2_2 The second element of mean for Mixture normal approximation of beta2.
#' @param prec2_2 The second element of precision vector for Mixture normal approximation of beta2.
#' @param N Number of replicates (default: 200).
#' @param n.burnin Number of burn-in samples for JAGS (default: 10000).
#' @param n.iter Total number of MCMC iterations for JAGS (default: 80000).
#' @param n.thin Thinning factor for JAGS (default: 20).
#' @param seed random seed.
#' 
#' @return A matrix containing MCMC samples of pediatric coefficients.

posterior_sample <- function(eps_H, lower, upper, n_ped, beta_adu, type,
                             weight_beta, weight_delta = NULL, beta_ped_all, 
                             w_1=seq(0.1,0.5,length=5), x_lower, x_upper, beta1_w, beta2_w, 
                             mu1_1, prec1_1, mu1_2, prec1_2, 
                             mu2_1, prec2_1, mu2_2, prec2_2, 
                             N = 200, n.burnin = 10000, n.iter = 80000, n.thin = 20,seed=1) {
  set.seed(seed)
  model.const = function(){
    for (i in 1:n_ped) {
      Y[i]~dbin(ilogit(coef_ped[1]+coef_ped[2]*X[i]), 1)
    }
    ind1~dcat(p1)
    latent1[1]~dnorm(0,0.001)
    latent1[2]~dnorm(mu1_1,prec1_1)
    latent1[3]~dnorm(mu1_2,prec1_2)
    coef_ped[1]<-latent1[ind1]
    
    ind2~dcat(p2)
    latent2[1]~dnorm(0,0.001)
    latent2[2]~dnorm(mu2_1,prec2_1)
    latent2[3]~dnorm(mu2_2,prec2_2)
    coef_ped[2]<-latent2[ind2]
  }
  
  # Validate input arguments
  if (type == "type I" && length(weight_beta) != 10) {
    stop("For type I, weight_beta must have exactly 10 elements.")
  }
  if (type == "power" && length(weight_beta) != 11) {
    stop("For power, weight_beta must have exactly 11 elements.")
  }
  if (type == "power" && is.null(weight_delta)) {
    stop("For power, weight_delta must be provided.")
  }
  if (type == "type I" && !"matrix" %in% class(beta_ped_all)) {
    stop("For type I, beta_ped_all must be a matrix.")
  }
  
  # Define MCMC sample storage
  n.row <- (n.iter - n.burnin) / n.thin
  temp <- matrix(0, nrow = N * n.row, ncol = 2)
  MCMC <- matrix(0, nrow = length(w_1) * N * n.row, ncol = 2)
  
  set.seed(1)  # Set seed for reproducibility
  
  # Iterate over w_1 values
  for (j in 1:length(w_1)) {
    p1 <- c(1 - w_1[j], w_1[j] * beta1_w)
    p2 <- c(1 - w_1[j], w_1[j] * beta2_w)
    
    for (i in 1:N) {
      print(paste0("i=", i, ", j=", j))
      
      
        repeat {
          repeat {
            if (type == "type I") {
              # Sample from beta_ped_all using weight_beta
              s <- sample(seq(1, 10), 1, replace = TRUE, prob = weight_beta)
              beta_ped <- as.vector(beta_ped_all[s,])
              
            } else if (type == "power") {
              # Sample delta index first
              del <- sample(seq(1, 10), 1, replace = TRUE, prob = weight_delta)
              temp_beta <- beta_ped_all[[del]]
              
              # Sample beta index from temp_beta
              s <- sample(seq(1, 11), 1, replace = TRUE, prob = weight_beta)
              beta_ped <- as.vector(temp_beta[s, ])
            }
            
            # Generate random x values within the given range
            x_ped <- runif(n_ped, x_lower, x_upper)
            
            # Compute p_gen and ensure values are within [0,1]
            q <- beta_ped[1] + x_ped * beta_ped[2]
            p_gen <- exp(q) / (1 + exp(q))
            p_gen <- pmax(0, pmin(1, p_gen))  # Constrain between [0,1]
            
            # Generate response values
            ped <- rbinom(n_ped, 1, prob = p_gen)
            
            # Create data frame
            data_ped <- data.frame(x = x_ped, y = ped)
            
            # Fit logistic regression model
            freq <- try(glm(y ~ x, data = data_ped, family = "binomial"), silent = TRUE)
            
            if (!inherits(freq, "try-error")) break
          }
          
          # Prepare MCMC data
          data_final <- list(Y = data_ped$y, X = data_ped$x, n_ped = n_ped, 
                             p1 = p1, p2 = p2, 
                             mu1_1 = mu1_1, prec1_1 = prec1_1, 
                             mu1_2 = mu1_2, prec1_2 = prec1_2,
                             mu2_1 = mu2_1, prec2_1 = prec2_1,
                             mu2_2 = mu2_2, prec2_2 = prec2_2)
          
          inits <- list(list("latent1" = rep(freq$coefficients[1], 3),
                             "latent2" = rep(freq$coefficients[2], 3)))
          
          # Run MCMC with JAGS
          AP <- try(jags(data = data_final, inits = inits, model = model.const,
                         n.chain = 1, n.burnin = n.burnin, n.thin = n.thin, 
                         parameters.to.save = c("coef_ped"), n.iter = n.iter),
                    silent = TRUE)
          
          if (!inherits(AP, "try-error")) break
        }
        
      
      # Store MCMC results
      temp[((i - 1) * n.row + 1):(i * n.row), 1:2] <- as.matrix(as.mcmc(AP)[[1]][, 1:2])
    }
    
    # Store all MCMC results
    MCMC[((j - 1) * N * n.row + 1):(j * N * n.row), 1:2] <- temp
  }
  
  return(MCMC)
}

###Example
#power

MCMC_results_power <- posterior_sample(
  eps_H = 0.2, lower = 2.5, upper = 5, 
  n_ped = 40, beta_adu = c(-2.830380, 1.412872),
  type = "power", weight_beta = c(1,1,1,1,1,1,2,2,2,4,8),
  weight_delta = c(1,2,4,8,4,2,1,1,1,1),
  beta_ped_all = result_power_final, w_1 = seq(0.1,0.5, length = 5),
  x_lower = 0, x_upper = 5,
  beta1_w = result$beta1_w, beta2_w = result$beta2_w, N=1,
  mu1_1 =result$mu1_1 , prec1_1 = result$prec1_1, mu1_2 = result$mu1_2, prec1_2 = result$prec1_2, 
  mu2_1 = result$mu2_1, prec2_1 = result$prec2_1, mu2_2 = result$mu2_2, prec2_2 = result$prec2_2
)

#type I

MCMC_results_typeI <- posterior_sample(
  eps_H = 0.2, lower = 2.5, upper = 5, 
  n_ped = 40, beta_adu = c(-2.830380, 1.412872),
  type = "type I", weight_beta = c(8,4,2,2,2,1,1,1,1,1),
  beta_ped_all = result_typeI, w_1 = seq(0.1,0.5, length = 5),
  x_lower = 0, x_upper = 5,
  beta1_w = result$beta1_w, beta2_w = result$beta2_w, N=1,
  mu1_1 =result$mu1_1 , prec1_1 = result$prec1_1, mu1_2 = result$mu1_2, prec1_2 = result$prec1_2, 
  mu2_1 = result$mu2_1, prec2_1 = result$prec2_1, mu2_2 = result$mu2_2, prec2_2 = result$prec2_2
)
