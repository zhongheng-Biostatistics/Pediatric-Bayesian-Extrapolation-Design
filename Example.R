source("typeI_coeff.r")
source("power_coeff.r")
source("Repp.r")
source("posterior_sample.r")
source("output.r")

eps_H = 0.2
lower = 2.5
upper = 5
n_adu=350
n_ped = 40
beta_adu = c(-2.830380, 1.412872)

## when H_0 holds, we obtain all the coefficients with respect to differnt eta
result_typeI <- typeI_coeff(eps_H,n_adu,lower,upper,beta_adu=beta_adu)
result_typeI<-result_typeI[[1]]

## when H_1 holds, we obtain all the coefficients with respect to differnt eta combined with delta
result_power <- power_coeff(eps_H,n_adu,lower,upper,beta_adu=beta_adu)
result_power_final<-list(result_power[[1]],result_power[[2]],result_power[[3]],result_power[[4]],result_power[[5]],result_power[[6]],result_power[[7]],result_power[[8]],result_power[[9]],result_power[[10]])  
## result_power_final combine all the pedatric coefficients sorted by the value of delta


###The prior of the pediatric coefficients
x <- c(lower, lower + 0.25 * (upper - lower), lower + 0.75 * (upper - lower))  
bias <- rep(0, length(x))
sd <- rep(0.01, length(x))

result <- Repp(beta_adu=beta_adu, lower=lower, upper=upper, x=x, bias=bias, sd=sd)


## when H_0 holds
weight_beta = c(8,4,2,2,2,1,1,1,1,1)  # the weight of different beta with respect to eta
w_1 = seq(0.1,0.5, length = 5)        # the weight of adult information in REPP
x_lower=0                             # Lower bound of the exposure range                       
x_upper=5                             # Upper bound of the exposure range
N=1                                   # Number of replicates


MCMC_results_typeI <- posterior_sample(
  eps_H = eps_H, lower = lower, upper = upper, 
  n_ped =n_ped, beta_adu = beta_adu,
  type = "type I", weight_beta = weight_beta,
  beta_ped_all = result_typeI, w_1 = w_1,
  x_lower = x_lower, x_upper = x_upper,
  beta1_w = result$beta1_w, beta2_w = result$beta2_w, N=N,
  mu1_1 =result$mu1_1 , prec1_1 = result$prec1_1, mu1_2 = result$mu1_2, prec1_2 = result$prec1_2, 
  mu2_1 = result$mu2_1, prec2_1 = result$prec2_1, mu2_2 = result$mu2_2, prec2_2 = result$prec2_2
)


## When H_1 holds
weight_beta = c(1,1,1,1,1,1,2,2,2,4,8)
weight_delta = c(1,2,4,8,4,2,1,1,1,1)

MCMC_results_power <- posterior_sample(
  eps_H = eps_H, lower = lower, upper = upper, 
  n_ped = n_ped, beta_adu = beta_adu,
  type = "power", weight_beta = weight_beta,
  weight_delta = weight_delta,
  beta_ped_all = result_power_final, w_1 = w_1,
  x_lower = x_lower, x_upper = x_upper,
  beta1_w = result$beta1_w, beta2_w = result$beta2_w, N=N,
  mu1_1 =result$mu1_1 , prec1_1 = result$prec1_1, mu1_2 = result$mu1_2, prec1_2 = result$prec1_2, 
  mu2_1 = result$mu2_1, prec2_1 = result$prec2_1, mu2_2 = result$mu2_2, prec2_2 = result$prec2_2
)


## Final output

post_eps = c(0.8, 0.85, 0.9, 0.95, 0.99)

# when H_0 holds
output_results_typeI <- output(eps_H=eps_H,
                         lower = lower, upper = upper, w_1 = w_1, 
                         N =N, beta_adu = beta_adu,
                         post_eps = post_eps, 
                         posterior_sample = MCMC_results_typeI,
                         n.burnin = 10000, n.iter = 80000, n.thin = 20,
                         type = "type I")

# when H_1 holds

output_results_power <- output(eps_H=eps_H,
                         lower = lower, upper = upper, w_1 = w_1, 
                         N =N, beta_adu = beta_adu,
                         post_eps = post_eps, 
                         posterior_sample = MCMC_results_power,
                         n.burnin = 10000, n.iter = 80000, n.thin = 20,
                         type = "power")