source("Simulation_function.R")


### Cross-sectional
cross_sectional_MC = function(n, m, M, G) {
  mu = rep(0, m)
  
  formula = h ~ 1
  
  beta1 = matrix(nrow = M, ncol = m)
  beta2 = matrix(nrow = M, ncol = m)
  vbeta = array(dim = c(M, 2*m, 2*m))
  
  set.seed(1)
  for (MC in 1:M) {
    data = subsample(G = G, n = n)
    
    h1.sum = sum(data$h1)
    h2.sum = sum(data$h2)
    beta10 = qlogis(h1.sum/h2.sum)
    beta20 = qlogis(h2.sum/nrow(data))
    
    fit.ugee = ugee_triangle2(formula = formula,
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

n = 250
m = 1
M = 1000

## transitivity = 0.4
adjustment_level = 0.002
true_p = true_p_triangle(sigma, xi, delta, type = "cs", adjustment_level = adjustment_level)

(true_trans = true_p$numerator / true_p$denominator)

beta10 = qlogis(true_trans)
beta20 = qlogis(true_p$denominator)

cs_n250_adj0_002_result = cross_sectional_MC(n, m, M, true_p$g)
save(cs_n250_adj0_002_result, file = "cross_sectional_n250_M1000_adj0_002.RData")

## transitivity = 0.2
adjustment_level = 0.008
true_p = true_p_triangle(sigma, xi, delta, type = "cs", adjustment_level = adjustment_level)

(true_trans = true_p$numerator / true_p$denominator)

beta10 = qlogis(true_trans)
beta20 = qlogis(true_p$denominator)

cs_n250_adj0_008_result = cross_sectional_MC(n, m, M, true_p$g)
save(cs_n250_adj0_008_result, file = "cross_sectional_n250_M1000_adj0_008.RData")

