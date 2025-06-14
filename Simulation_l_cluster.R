source("Simulation_function.R")


### longitudinal
longitudinal_MC = function(n, m, M, sigma, xi, delta, beta10, beta20, pi0) {
  mu = rep(0, m)
  
  formula = w * h ~ 1
  
  beta1 = matrix(nrow = M, ncol = m)
  beta2 = matrix(nrow = M, ncol = m)
  pi = matrix(nrow = M, ncol = m-1)
  vtau = array(dim = c(M, 3*m-1, 3*m-1))
  
  set.seed(1)
  for (MC in 1:M) {
    restart = 1
    while(restart == 1) {
      data = network_data(n = n, m = m,
                          mu = mu, sigma = sigma, xi = xi,
                          delta = delta,
                          pi0 = pi0)
      
      trans1 = sum(floor(data$w_1) * data$h1) / sum(floor(data$w_1) * data$h2)
      trans2 = sum(floor(data$w_2) * data$h1) / sum(floor(data$w_2) * data$h2)
      trans3 = sum(data$h1) / sum(data$h2)
      
      if (!any(c(trans1, trans2, trans3) %in% c(0, 1))) {
        restart = 0
        fit.ugee = ugee_triangle(formula = formula,
                                 data = data,
                                 n = n,
                                 id = c("i", "j", "k"),
                                 corstr = "independence",
                                 maxit = 200,
                                 tol = 0.001,
                                 true_param = c(beta10, beta20, pi0))
      }
    }
    
    beta1[MC, ] = fit.ugee$beta1
    beta2[MC, ] = fit.ugee$beta2
    pi[MC, ] = fit.ugee$pi
    vtau[MC, , ] = fit.ugee$vtau
  }
  
  return(list(beta1 = beta1,
              beta2 = beta2,
              pi = pi,
              vtau = vtau))
}


# Simulated transitivity
sigma = c(12, 3, 2)
xi = c(-3, -1, -0.5)
delta = 1.5
pi0 = c(0.4, 0.7, 1)

true_p = true_p_triangle(sigma, xi, delta, pi0, type = "l")

(true_trans = true_p$numerator / true_p$denominator)

beta10 = qlogis(true_trans)
beta20 = qlogis(true_p$denominator)

## n = 150, m = 3, M = 1000
n = 150
m = 3
M = 1000

l_m3_n150_result = longitudinal_MC(n, m, M, sigma, xi, delta, beta10, beta20, pi0)
save(l_m3_n150_result, file = "longitudinal_m3_n150_M1000.RData")

## n = 250, m = 3, M = 1000
n = 250
l_m3_n250_result = longitudinal_MC(n, m, M, sigma, xi, delta, beta10, beta20, pi0)
save(l_m3_n250_result, file = "longitudinal_m3_n250_M1000.RData")

## n = 500, m = 3, M = 1000
n = 500
l_m3_n500_result = longitudinal_MC(n, m, M, sigma, xi, delta, beta10, beta20, pi0)
save(l_m3_n500_result, file = "longitudinal_m3_n500_M1000.RData")

