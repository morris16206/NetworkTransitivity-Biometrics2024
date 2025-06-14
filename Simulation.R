rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


source("Simulation_function.R")

library(qqplotr)
library(gtable)
library(gridExtra)

### Parameters for simulated transitivity
real.sigma = c(12, 3, 2)
real.xi = c(-3, -1, -0.5)

### Cross-sectional
## n = 150/250/500, m = 1, M = 1000
n = 150
m = 1
M = 1000

# Simulated transitivity
sigma = 2
xi = -0.5
delta = 1.5

true_p = true_p_triangle(sigma, xi, delta, type = "cs")

(true_trans = true_p$numerator / true_p$denominator)

beta10 = qlogis(true_trans)

load("cross_sectional_n150_M1000.RData")
#load("cross_sectional_n250_M1000.RData")
#load("cross_sectional_n500_M1000.RData")

cs_result = cs_n150_result
#cs_result = cs_n250_result
#cs_result = cs_n500_result

# Estimated f(beta)
fbeta = matrix(nrow = M, ncol = m)
for (MC in 1:M) {
  fbeta[MC, ] = expit(cs_result$beta1[MC, ])
}
colMeans(fbeta)

# Asymptotic variance of f(beta)
v.fbeta = matrix(nrow = M, ncol = m)
for (MC in 1:M) {
  diff = exp(cs_result$beta1[MC, ]) / (1 + exp(cs_result$beta1[MC, ]))^2
  v.fbeta[MC, ] = cs_result$vbeta[MC, 1, 1] * diff^2
}
colMeans(v.fbeta) / n

# Empirical variance of f(beta)
apply(fbeta, 2, var)

# Empirical type I error
alpha.Wald = c()
for (MC in 1:M) {
  Wald = t(cs_result$beta1[MC, ] - beta10) %*% ginv(cs_result$vbeta[MC, 1, 1] / n) %*% (cs_result$beta1[MC, ] - beta10)
  alpha.Wald = c(alpha.Wald, Wald >= qchisq(0.95, 1))
}
mean(alpha.Wald)

# qq-plot
qq_cs_n150 = ggplot(data.frame(fbeta), aes(sample = fbeta)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 150)") +
  scale_x_continuous(breaks = seq(0.74, 0.88, 0.02)) +
  scale_y_continuous(breaks = seq(0.74, 0.88, 0.02)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_cs_n250 = ggplot(data.frame(fbeta), aes(sample = fbeta)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 250)") +
  scale_x_continuous(breaks = seq(0.76, 0.86, 0.02)) +
  scale_y_continuous(breaks = seq(0.76, 0.86, 0.02)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_cs_n500 = ggplot(data.frame(fbeta), aes(sample = fbeta)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 500)") +
  scale_x_continuous(breaks = seq(0.76, 0.86, 0.02)) +
  scale_y_continuous(breaks = seq(0.76, 0.86, 0.02)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

grid.arrange(
  qq_cs_n150,
  qq_cs_n250,
  qq_cs_n500,
  ncol = 3
)

## transitivity = 0.4
adjustment_level = 0.002
true_p = true_p_triangle(sigma, xi, delta, type = "cs", adjustment_level = adjustment_level)

(true_trans = true_p$numerator / true_p$denominator)

beta10 = qlogis(true_trans)

load("cross_sectional_n250_M1000_adj0_002.RData")

cs_result = cs_n250_adj0_002_result

# Estimated f(beta)
fbeta = matrix(nrow = M, ncol = m)
for (MC in 1:M) {
  fbeta[MC, ] = expit(cs_result$beta1[MC, ])
}
colMeans(fbeta)

# Empirical variance of f(beta)
apply(fbeta, 2, var) / (1 - n / 1400)

## transitivity = 0.2
adjustment_level = 0.008
true_p = true_p_triangle(sigma, xi, delta, type = "cs", adjustment_level = adjustment_level)

(true_trans = true_p$numerator / true_p$denominator)

beta10 = qlogis(true_trans)

load("cross_sectional_n250_M1000_adj0_008.RData")

cs_result = cs_n250_adj0_008_result

# Estimated f(beta)
fbeta = matrix(nrow = M, ncol = m)
for (MC in 1:M) {
  fbeta[MC, ] = expit(cs_result$beta1[MC, ])
}
colMeans(fbeta)

# Empirical variance of f(beta)
apply(fbeta, 2, var) / (1 - n / 1400)


### Longitudinal
## n = 150/250/500, m = 3, M = 1000
n = 150
m = 3
M = 1000

# Simulated transitivity
sigma = real.sigma
xi = real.xi
delta = 1.5
nonmissing = c(0.4, 0.7, 1)

true_p = true_p_triangle(sigma = sigma, xi = xi, delta = delta, pi0 = nonmissing, type = "l")

(true_trans = true_p$numerator / true_p$denominator)

beta10 = qlogis(true_trans)

load("longitudinal_m3_n150_M1000.RData")
#load("longitudinal_m3_n250_M1000.RData")
#load("longitudinal_m3_n500_M1000.RData")

l_result = l_m3_n150_result
#l_result = l_m3_n250_result
#l_result = l_m3_n500_result

# Estimated f(beta)
fbeta = matrix(nrow = M, ncol = m)
for (MC in 1:M) {
  fbeta[MC, ] = expit(l_result$beta1[MC, ])
}
colMeans(fbeta)

# Asymptotic variance of f(beta)
v.fbeta = matrix(nrow = M, ncol = m)
#v.fbeta = array(dim = c(M, m, m))
for (MC in 1:M) {
  grad = exp(l_result$beta1[MC, ]) / (1 + exp(l_result$beta1[MC, ]))^2
  v.fbeta[MC, ] = diag(l_result$vtau[MC, , ])[1:m] * grad^2
  
  #v.fbeta[MC, , ] = diag(grad) %*% l_result$vbeta[MC, , ] %*% diag(grad)
  #v.fbeta[MC, ] = diag(l_result$vtau[MC, 1:m, 1:m]) * (grad * nonmissing^3)^2
}
colMeans(v.fbeta) / n
#apply(v.fbeta, 2:3, mean) / n

# Empirical variance of f(beta)
apply(fbeta, 2, var)
#fbeta.center = sweep(fbeta, 2, colMeans(fbeta))
#t(fbeta.center) %*% fbeta.center / (M-1)

# Empirical type I errors
alpha.Wald = matrix(nrow = M, ncol = m)
alpha.Wald.omnibus = matrix(nrow = M, ncol = 1)
L = matrix(0, nrow = 2, ncol = 3)
L[, 1] = -1
L[row(L) == col(L) - 1] = 1
for (MC in 1:M) {
  Wald = t(l_result$beta1[MC, ] - beta10) * (n / diag(l_result$vtau[MC, , ])[1:m]) * (l_result$beta1[MC, ] - beta10)
  alpha.Wald[MC, ] = (Wald >= qchisq(0.95, 1))
  
  Wald.omnibus = n * t(L %*% l_result$beta1[MC, ]) %*% solve(L %*% l_result$vtau[MC, 1:m, 1:m] %*% t(L)) %*% (L %*% l_result$beta1[MC, ])
  alpha.Wald.omnibus[MC, ] = (Wald.omnibus >= qchisq(0.95, 2))
}
colMeans(alpha.Wald)
colMeans(alpha.Wald.omnibus)

# qq-plot
qq_l_n150_1 = ggplot(data.frame(fbeta), aes(sample = X1)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 60, time = 1)") +
  scale_x_continuous(breaks = seq(0.68, 1, 0.08)) +
  scale_y_continuous(breaks = seq(0.68, 1, 0.08)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_l_n150_2 = ggplot(data.frame(fbeta), aes(sample = X2)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 105, time = 2)") +
  scale_x_continuous(breaks = seq(0.72, 1, 0.04)) +
  scale_y_continuous(breaks = seq(0.72, 1, 0.04)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_l_n150_3 = ggplot(data.frame(fbeta), aes(sample = X3)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 150, time = 3)") +
  scale_x_continuous(breaks = seq(0.72, 0.83, 0.04)) +
  scale_y_continuous(breaks = seq(0.72, 0.84, 0.04)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_l_n250_1 = ggplot(data.frame(fbeta), aes(sample = X1)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 100, time = 1)") +
  scale_x_continuous(breaks = seq(0.76, 1, 0.08)) +
  scale_y_continuous(breaks = seq(0.76, 1, 0.08)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_l_n250_2 = ggplot(data.frame(fbeta), aes(sample = X2)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 175, time = 2)") +
  scale_x_continuous(breaks = seq(0.72, 1, 0.04)) +
  scale_y_continuous(breaks = seq(0.72, 1, 0.04)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_l_n250_3 = ggplot(data.frame(fbeta), aes(sample = X3)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 250, time = 3)") +
  scale_x_continuous(breaks = seq(0.72, 1, 0.04)) +
  scale_y_continuous(breaks = seq(0.72, 1, 0.04)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_l_n500_1 = ggplot(data.frame(fbeta), aes(sample = X1)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 200, time = 1)") +
  scale_x_continuous(breaks = seq(0.76, 1, 0.08)) +
  scale_y_continuous(breaks = seq(0.76, 1, 0.08)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_l_n500_2 = ggplot(data.frame(fbeta), aes(sample = X2)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 350, time = 2)") +
  scale_x_continuous(breaks = seq(0.72, 1, 0.04)) +
  scale_y_continuous(breaks = seq(0.72, 1, 0.04)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

qq_l_n500_3 = ggplot(data.frame(fbeta), aes(sample = X3)) +
  stat_qq_point(distribution = "norm") +
  stat_qq_line(distribution = "norm", color = "blue") +
  stat_qq_band(distribution = "norm",
               bandType = "pointwise",
               alpha = 0,
               linetype = "dashed",
               color = "blue") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "QQ-Plot (n = 500, time = 3)") +
  scale_x_continuous(breaks = seq(0.72, 1, 0.04)) +
  scale_y_continuous(breaks = seq(0.72, 1, 0.04)) +
  theme_bw() +
  theme(title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

grid.arrange(
  arrangeGrob(qq_l_n150_1, qq_l_n150_2, qq_l_n150_3, nrow = 1),
  arrangeGrob(qq_l_n250_1, qq_l_n250_2, qq_l_n250_3, nrow = 1),
  arrangeGrob(qq_l_n500_1, qq_l_n500_2, qq_l_n500_3, nrow = 1),
  nrow = 3
)

