#' Return the pediatric coefficient when H_1 holds 
#' 
#' This function will give corresponding pediatric coefficient for eta=0.1,0.15,0.25,...,0.95 combine with delta=eps_H*k/10,k=0,1...9.
#' 
#' @param eps_H Threshold in H_0.
#' @param n_adu sample size of adult dataset in the literature.
#' @param lower Lower bound of the region of interest.
#' @param upper Upper bound of the region of interest.
#' @param beta_adu Adult's coefficients.

#' @return A matrix containing all the pediatric coefficients with respect to eta=0.1,0.15,0.25,...,0.95 combine with delta=eps_H*k/10,k=0,1...9.


library(nleqslv)
library(dplyr)

power_coeff <- function(eps_H, n_adu, lower, upper, beta_adu) {
  set.seed(10)
  eps_1_seq <- eps_H / 10 * seq(0, 9)  
  coefficients_list <- list()  
  
  fit <- function(x, beta) {
    exp(beta[1] + x * beta[2]) / (1 + exp(beta[1] + x * beta[2]))
  }
  
  max_dev <- function(x, beta1, beta2) {
    abs(fit(x, beta1) - fit(x, beta2))
  }
  
  eq_1 <- function(x, para1 = c(lower, upper), y) {
    c(fit(para1[1], x) - y[1], fit(para1[2], x) - y[2])
  }
  
  for (eps_1 in eps_1_seq) {
    if (eps_1 == 0) {
      beta_matrix <- matrix(rep(beta_adu, 11), nrow = 11, byrow = TRUE)
      coefficients_list[[as.character(eps_1)]] <- beta_matrix
      next
    }
    
    y_resp <- c(max(fit(lower, beta_adu) - eps_1, 0), fit(upper, beta_adu))
    sol_1 <- nleqslv(beta_adu, eq_1, y = y_resp)
    
    y_resp_2 <- c(fit(lower, beta_adu), fit(upper, beta_adu) - eps_1)
    sol_2 <- nleqslv(beta_adu, eq_1, y = y_resp_2)
    
    lower_mark <- sol_1$x[1] + sol_1$x[2] * lower
    upper_mark <- sol_2$x[1] + sol_2$x[2] * upper
    
    n <- 500  
    t1 <- seq(lower_mark, beta_adu[1] + beta_adu[2] * lower, length = n)
    t2 <- seq(upper_mark, beta_adu[1] + beta_adu[2] * upper, length = n)
    
    result <- matrix(0, ncol = 4, nrow = n)
    colnames(result) <- c("difference_area", "max_difference", "beta_0", "beta_1")
    
    beta_temp <- c(0, 0)
    
    for (j in 1:n) {
      temp <- numeric(n)
      for (i in 1:n) {
        beta_temp[2] <- (t1[i] - t2[j]) / (lower - upper)
        beta_temp[1] <- t1[i] - beta_temp[2] * lower
        temp[i] <- max(
          optimize(max_dev, c(lower, upper), maximum = TRUE, beta1 = beta_adu, beta2 = beta_temp)$objective,
          max_dev(lower, beta_adu, beta_temp),
          max_dev(upper, beta_adu, beta_temp)
        )
      }
      ind <- which.min(abs(temp - eps_1))
      result[j, ] <- c(
        integrate(max_dev, lower, upper, beta1 = beta_adu, beta2 = beta_temp)$value,
        temp[ind], 
        t1[ind] - (t1[ind] - t2[j]) / (lower - upper) * lower,
        (t1[ind] - t2[j]) / (lower - upper)
      )
    }
    
    tol <- 1e-4
    result <- as.data.frame(result)
    result_final <- filter(result, abs(max_difference - eps_1) < tol)
    
    result_final$ratio <- result_final[, 1] / max(result_final[, 1])
    
    ind_10 <- which.min(abs(result_final$ratio - 0.1))
    ind_15 <- which.min(abs(result_final$ratio - 0.15))
    ind_25 <- which.min(abs(result_final$ratio - 0.25))
    ind_35 <- which.min(abs(result_final$ratio - 0.35))
    ind_45 <- which.min(abs(result_final$ratio - 0.45))
    ind_55 <- which.min(abs(result_final$ratio - 0.55))
    ind_65 <- which.min(abs(result_final$ratio - 0.65))
    ind_75 <- which.min(abs(result_final$ratio - 0.75))
    ind_85 <- which.min(abs(result_final$ratio - 0.85))
    ind_95 <- which.min(abs(result_final$ratio - 0.95))
    ind_max <- which.max(result_final[, 1])
    
    beta_10 <- as.numeric(result_final[ind_10, c(3, 4)])
    beta_15 <- as.numeric(result_final[ind_15, c(3, 4)])
    beta_25 <- as.numeric(result_final[ind_25, c(3, 4)])
    beta_35 <- as.numeric(result_final[ind_35, c(3, 4)])
    beta_45 <- as.numeric(result_final[ind_45, c(3, 4)])
    beta_55 <- as.numeric(result_final[ind_55, c(3, 4)])
    beta_65 <- as.numeric(result_final[ind_65, c(3, 4)])
    beta_75 <- as.numeric(result_final[ind_75, c(3, 4)])
    beta_85 <- as.numeric(result_final[ind_85, c(3, 4)])
    beta_95 <- as.numeric(result_final[ind_95, c(3, 4)])
    beta_max <- as.numeric(result_final[ind_max, c(3, 4)])
    
    beta_matrix <- rbind(beta_10, beta_15, beta_25, beta_35, beta_45, beta_55,
                         beta_65, beta_75, beta_85, beta_95, beta_max)
    
    coefficients_list[[as.character(eps_1)]] <- beta_matrix
  }
  
  return(coefficients_list)
}

eps_H <- 0.2
n_adu=350
lower=2.5
upper=5

# Example
Time1<-Sys.time()
result_power <- power_coeff(eps_H,n_adu,lower,upper,beta_adu=c(-2.830380,1.412872))
Sys.time()-Time1

# eps_1=0.04 with its coefficient
result_power_final<-list(result_power[[1]],result_power[[2]],result_power[[3]],result_power[[4]],result_power[[5]],result_power[[6]],result_power[[7]],result_power[[8]],result_power[[9]],result_power[[10]])  
