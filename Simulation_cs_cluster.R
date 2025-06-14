source("Simulation_function.R")


### Cross-sectional
cross_sectional_MC = function(n, m, M, sigma, xi, delta, beta10, beta20) {
  mu = rep(0, m)
  
  formula = h ~ 1
  
  beta1 = matrix(nrow = M, ncol = m)
  beta2 = matrix(nrow = M, ncol = m)
  vbeta = array(dim = c(M, 2*m, 2*m))
  
  set.seed(1)
  for (MC in 1:M) {
    data = network_data(n = n, m = m,
                        mu = mu, sigma = sigma, xi = xi,
                        delta = delta)
    
    fit.ugee = ugee_triangle(formula = formula,
                             data = data,
                             n = n,
                             id = c("i", "j", "k"),
                             corstr = "independence",
                             maxit = 200,
                             tol = 0.001,
                             true_param = c(beta10, beta20))
    
    beta1[MC, ] = fit.ugee$beta1
    beta2[MC, ] = fit.ugee$beta2
    vbeta[MC, , ] = fit.ugee$vbeta
  }
  
  return(list(beta1 = beta1,
              beta2 = beta2,
              vbeta = vbeta))
}


# Simulated transitivity
sigma = 2
xi = -0.5
delta = 1.5

true_p = true_p_triangle(sigma, xi, delta, type = "cs")

(true_trans = true_p$numerator / true_p$denominator)

beta10 = qlogis(true_trans)
beta20 = qlogis(true_p$denominator)

## n = 150, m = 1, M = 1000
n = 150
m = 1
M = 1000

cs_n150_result = cross_sectional_MC(n, m, M, sigma, xi, delta, beta10, beta20)
save(cs_n150_result, file = "cross_sectional_n150_M1000.RData")

## n = 250, m = 1, M = 1000
n = 250
cs_n250_result = cross_sectional_MC(n, m, M, sigma, xi, delta, beta10, beta20)
save(cs_n250_result, file = "cross_sectional_n250_M1000.RData")

## n = 500, m = 1, M = 1000
n = 500
cs_n500_result = cross_sectional_MC(n, m, M, sigma, xi, delta, beta10, beta20)
save(cs_n500_result, file = "cross_sectional_n500_M1000.RData")

