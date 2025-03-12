library(R2jags)
library(R2WinBUGS)
library(flexmix)

#' @param beta_adu: Adult's coefficients
#' @param lower: The lower bound of the interested interval
#' @param upper: The upper bound of the interested interval
#' @param x: Selected points which have the response information (vector of length 2, 3, or 4)
#' @param bias: The initial guess of the difference of response between pediatric and adult at each x (same length as x)
#' @param sd: The standard deviation of the bias at each x (same length as x)
#' @param n.replicate: Number of replications for each x
#' @param n.burnin: Number of burn-in iterations in JAGS
#' @param n.iter: Total number of iterations in JAGS
#' @param n.thin: Thinning interval in JAGS
#' @param inits: Initial values of beta for JAGS

Repp <- function(beta_adu, lower, upper, x, bias, sd,
                 n.replicate = 1000, n.burnin = 10000, n.iter = 60000,
                 n.thin = 20, inits = list(list('coef_h1' = c(-1, 1)))) {

  # Get the number of selected points (x dimension)
  x_dim <- length(x)

  # Generate simulated data
  set.seed(123)
  y.sim <- matrix(0, ncol = x_dim, nrow = n.replicate)

  for (i in 1:n.replicate) {
    o <- beta_adu[1] + beta_adu[2] * x+ rnorm(x_dim, bias, sd)
    p_gen=exp(o)/(1+exp(o))
    y <- rbinom(x_dim, 1, prob = p_gen)
    y.sim[i, ] <- y
  }

  # Create the data matrix with responses and corresponding x values
  data_H <- do.call(cbind, lapply(1:x_dim, function(i) cbind(y.sim[, i], rep(x[i], n.replicate))))

  # Define the BUGS model dynamically based on the x dimension
  model_H <- function() {
    for (i in 1:n.replicate) {
      for (j in 1:x_dim) {
        temp[i, j] <- ilogit(coef_h1[1] + X[i, j] * coef_h1[2])
      }
    }
    for (i in 1:n.replicate) {
      for (j in 1:x_dim) {
        Y[i, j] ~ dbin(temp[i, j], 1)
      }
    }
    coef_h1[1] ~ dnorm(0, 0.001)
    coef_h1[2] ~ dnorm(0, 0.001)
  }

  # Prepare the data list for JAGS
  data_final <- list(
    Y = as.matrix(data_H[, seq(1, 2 * x_dim, by = 2)]),  # Response matrix
    X = as.matrix(data_H[, seq(2, 2 * x_dim, by = 2)]),  # Corresponding x values
    n.replicate = n.replicate,
    x_dim = x_dim
  )

  # Run JAGS
  H_mcmc <- jags(
    data = data_final, inits = inits, model = model_H,
    n.chains = 1, n.burnin = n.burnin, n.thin = n.thin,
    parameters.to.save = c("coef_h1"), n.iter = n.iter
  )

  coef_H <- as.matrix(as.mcmc(H_mcmc))[, 1:2]
  colnames(coef_H) <- c("coef1", "coef2")

  # Use a mixture model to approximate the posterior
  mo1 <- FLXMRglm(family = "gaussian")
  mo2 <- FLXMRglm(family = "gaussian")

  flexfit_beta1 <- flexmix(coef1 ~ 1, data = as.data.frame(coef_H), k = 2, model = list(mo1, mo2))
  beta1_w<-table(clusters(flexfit))
  beta1_w<-beta1_w/sum(beta1_w)
  beta1_1 <- parameters(flexfit_beta1, component = 1)[[1]]
  beta1_2 <- parameters(flexfit_beta1, component = 2)[[1]]

  flexfit_beta2 <- flexmix(coef2 ~ 1, data = as.data.frame(coef_H), k = 2, model = list(mo1, mo2))
  beta2_w<-table(clusters(flexfit))
  beta2_w<-beta2_w/sum(beta2_w)
  beta2_1 <- parameters(flexfit_beta2, component = 1)[[1]]
  beta2_2 <- parameters(flexfit_beta2, component = 2)[[1]]

  # Compute beta1 parameters
  mu1_1 <- beta1_1[1]
  prec1_1 <- 1 / beta1_1[2]^2
  mu1_2 <- beta1_2[1]
  prec1_2 <- 1 / beta1_2[2]^2

  # Compute beta2 parameters
  mu2_1 <- beta2_1[1]
  prec2_1 <- 1 / beta2_1[2]^2
  mu2_2 <- beta2_2[1]
  prec2_2 <- 1 / beta2_2[2]^2

  return(list(
    mu1_1 = mu1_1, prec1_1 = prec1_1,
    mu1_2 = mu1_2, prec1_2 = prec1_2,
    beta1_w=beta1_w,
    mu2_1 = mu2_1, prec2_1 = prec2_1,
    mu2_2 = mu2_2, prec2_2 = prec2_2,
    beta2_w=beta2_w
  ))
}

# Example
beta_adu <- c(-2.830380, 1.412872)
lower <- 2.5
upper <- 5
x <- c(lower, lower + 0.25 * (upper - lower), lower + 0.75 * (upper - lower))
bias <- rep(0, length(x))
sd <- rep(0.01, length(x))

result <- Repp(beta_adu, lower, upper, x, bias, sd)

# Print results
print(result)


