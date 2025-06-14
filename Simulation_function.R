library(EnvStats)
library(igraph)
library(RcppAlgos)
library(data.table)
library(plyr)
library(tidyverse)
library(locfit)


Rcpp::sourceCpp("edge_dist.cpp")
Rcpp::sourceCpp("alpha_estimation.cpp")
Rcpp::sourceCpp("Gauss_newton_triangle.cpp")
Rcpp::sourceCpp("asym_var_triangle.cpp")


### Function to approximate the true probability of forming a triangle
true_p_triangle = function(sigma, xi, delta, pi0 = rep(1, length(sigma)), type, adjustment_level = 0) {
  # Set seed for reproducibility
  set.seed(123)
  
  # Number of samples
  if (type == "cs") {
    num_samples <- 1400 # Cross-sectional
  } else if (type == "l") {
    num_samples <- 1500 # Longitudinal
  }
  
  # Create the proportion vector for each time point
  prop <- pi0 - lag(pi0)
  prop[1] <- pi0[1]
  
  # Vector to store the true probabilities for each time point
  true_p_numerator <- c()
  true_p_denominator <- c()
  
  # Loop over each time point
  for (time in 1:length(sigma)) {
    # Create a temporary vector to store the proportion up to the current time point
    prop.tmp <- prop.table(prop[1:time])
    
    # Generate samples from mixture generalized extreme value distributions
    samples <- c()
    for (time_iter in 1:time) {
      if (prop.tmp[time_iter] != 0) {
        samples <- c(samples, rgevd(round(num_samples * prop.tmp[time_iter]), 0, sigma[time_iter], xi[time_iter]))
      }
    }
    
    # Calculate the distance matrix
    dist_matrix <- as.matrix(dist(as.matrix(samples)))
    
    # Create the adjacency matrix
    adj_matrix <- (dist_matrix <= delta) + 0
    diag(adj_matrix) <- 0
    
    # Create the graph and reduce transitivity
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    
    # Rewire to decrease transitivity
    if (adjustment_level != 0) {
      g <- g %>%
        rewire(keeping_degseq(niter = 10^6))
      
      tri_vec <- triangles(g)
      tri_mat <- matrix(tri_vec, ncol = 3, byrow = T)
      n_triangles <- length(tri_vec) / 3
      n_adjusted <- floor(n_triangles * adjustment_level)
      
      tri_selected_ids <- sample(1:n_triangles, n_adjusted)
      
      edge_selected_mat <- matrix(NA_integer_, nrow = 2, ncol = n_adjusted)
      
      for (i in 1:n_adjusted) {
        tri = tri_mat[tri_selected_ids[i], ]
        edge_selected_mat[, i] = combn(tri, 2)[, sample(1:3, 1)]
      }
      
      edge_selected_ids <- get.edge.ids(g, edge_selected_mat)
      
      g <- delete_edges(g, edge_selected_ids)
      adj_matrix <- as_adjacency_matrix(g, sparse = F)
    }
    
    # Count triangles and potentially transitive triads
    tri_vec <- triangles(g)
    n_triangles <- length(tri_vec) / 3
    
    # Count potentially transitive triads
    adj_matrix2 <- adj_matrix %*% adj_matrix
    diag(adj_matrix2) <- 0
    n_potential_trans <- sum(adj_matrix2) / 2
    
    # Calculate the proportion of values less than or equal to delta and return it
    n_choose_3 <- comboCount(num_samples, 3)
    true_p_numerator <- c(true_p_numerator, n_triangles / n_choose_3)
    true_p_denominator <- c(true_p_denominator, n_potential_trans / n_choose_3 / 3)
  }
  
  # Return the vector of true probabilities for each time point
  return(list(numerator = true_p_numerator,
              denominator = true_p_denominator,
              g = g))
}

network_data = function(n, m,
                        mu, sigma, xi,
                        delta,
                        pi0 = 1) {
  # Sample observing indicator
  if (!identical(pi0, 1)) {
    prop = pi0 - lag(pi0)
    prop[1] = pi0[1]
    
    Z = rmultinom(n, 1, prob = prop)
    n_m = cumsum(rowSums(Z))
  }
  
  # Sample network data
  Y = matrix(nrow = n, ncol = 1)
  ord.index = c()
  if (m != 1) {
    for (time in 1:m) {
      Y.tmp = rgevd(n, location = mu[time], scale = sigma[time], shape = xi[time])
      Y = ifelse(Z[time, ] == 1, Y.tmp, Y)
      ord.index = c(ord.index, which(Z[time, ] == 1))
    }
    Y = as.matrix(Y[ord.index])
  } else {
    Y[1:n, ] = rgevd(n, location = mu, scale = sigma, shape = xi)
  }
  
  edge.index = comboGeneral(n, 2)
  
  edge.dist = edge_dist(edge.index, Y)
  
  h.pair = ifelse(edge.dist <= delta, 1, 0)
  
  networkData = setDT(as.data.frame(comboGeneral(n, 3)))
  
  if (m == 1) {
    g = graph_from_edgelist(matrix(edge.index[h.pair == 1, ], ncol = 2), directed = F)
    triangle.list = matrix(triangles(g), ncol = 3, byrow = T)
    
    if (identical(triangle.list, matrix(integer(0), ncol = 3))) {
      networkData[, "h1" := 0]
    } else {
      triangle.list = t(apply(triangle.list, 1, sort))
      triangle.list = triangle.list[do.call(order, lapply(1:ncol(triangle.list), function(i) triangle.list[, i])), ]
      triangle.list = cbind(matrix(triangle.list, ncol = 3), 1)
      adj.edge.list = lapply(get.adjlist(g), function(x) {
        if (length(x) > 1) {
          comboGeneral(as.vector(x), 2)
        }
      })
      length(adj.edge.list) = n
      adj.edge.list = lapply(1:n, function(i) {
        x = adj.edge.list[[i]]
        if (!is.null(x)) {
          x = cbind(x, i)
          t(apply(x, 1, sort))
        }
      })
      adj.edge.list = Filter(Negate(is.null), adj.edge.list)
      adj.edge.list = as.data.frame(do.call(rbind, adj.edge.list)) %>%
        group_by(across(everything())) %>%
        summarise(n = n() / 3, .groups = "drop")
      colnames(adj.edge.list) = paste0("V", 1:4)
      networkData[setDT(as.data.frame(triangle.list)), on = paste0("V", 1:3), "h1" := i.V4]
      networkData[setDT(adj.edge.list), on = paste0("V", 1:3), "h2" := i.V4]
      networkData[is.na(networkData[["h1"]]), "h1"] = 0
      networkData[is.na(networkData[["h2"]]), "h2"] = 0
    }
  } else {
    g.list = list()
    for (time in 1:m) {
      edge.index.time = edge.index[, 1] <= n_m[time] & edge.index[, 2] <= n_m[time]
      g = graph_from_edgelist(matrix(edge.index[edge.index.time == T & h.pair == 1, ], ncol = 2), directed = F)
      
      if (time < m) {
        triangle.index.time = ((networkData[, 1] <= n_m[time]) + (networkData[, 2] <= n_m[time]) + (networkData[, 3] <= n_m[time])) / 3
        networkData[, paste0("w_", time) := triangle.index.time]
      }
    }
    triangle.list = matrix(triangles(g), ncol = 3, byrow = T)
    if (identical(triangle.list, matrix(integer(0), ncol = 3))) {
      networkData[, "h1" := 0]
    } else {
      triangle.list = t(apply(triangle.list, 1, sort))
      triangle.list = triangle.list[do.call(order, lapply(1:ncol(triangle.list), function(i) triangle.list[, i])), ]
      triangle.list = cbind(matrix(triangle.list, ncol = 3), 1)
      adj.edge.list = lapply(get.adjlist(g), function(x) {
        if (length(x) > 1) {
          comboGeneral(as.vector(x), 2)
        }
      })
      length(adj.edge.list) = n
      adj.edge.list = lapply(1:n, function(i) {
        x = adj.edge.list[[i]]
        if (!is.null(x)) {
          x = cbind(x, i)
          t(apply(x, 1, sort))
        }
      })
      adj.edge.list = Filter(Negate(is.null), adj.edge.list)
      adj.edge.list = as.data.frame(do.call(rbind, adj.edge.list)) %>%
        group_by(across(everything())) %>%
        summarise(n = n() / 3, .groups = "drop")
      colnames(adj.edge.list) = paste0("V", 1:4)
      networkData[setDT(as.data.frame(triangle.list)), on = paste0("V", 1:3), "h1" := i.V4]
      networkData[setDT(adj.edge.list), on = paste0("V", 1:3), "h2" := i.V4]
      networkData[is.na(networkData[["h1"]]), "h1"] = 0
      networkData[is.na(networkData[["h2"]]), "h2"] = 0
    }
  }
  
  colnames(networkData)[1:3] = c("i", "j", "k")
  
  return(networkData)
}

subsample <- function(G, n) {
  sampled_nodes <- sample(V(G), size = n, replace = F)
  g <- induced_subgraph(G, sampled_nodes)
  
  networkData = setDT(as.data.frame(comboGeneral(n, 3)))
  
  triangle.list = matrix(triangles(g), ncol = 3, byrow = T)
  
  if (identical(triangle.list, matrix(integer(0), ncol = 3))) {
    networkData[, "h1" := 0]
  } else {
    triangle.list = t(apply(triangle.list, 1, sort))
    triangle.list = triangle.list[do.call(order, lapply(1:ncol(triangle.list), function(i) triangle.list[, i])), ]
    triangle.list = cbind(matrix(triangle.list, ncol = 3), 1)
    adj.edge.list = lapply(get.adjlist(g), function(x) {
      if (length(x) > 1) {
        comboGeneral(as.vector(x), 2)
      }
    })
    length(adj.edge.list) = n
    adj.edge.list = lapply(1:n, function(i) {
      x = adj.edge.list[[i]]
      if (!is.null(x)) {
        x = cbind(x, i)
        t(apply(x, 1, sort))
      }
    })
    adj.edge.list = Filter(Negate(is.null), adj.edge.list)
    adj.edge.list = as.data.frame(do.call(rbind, adj.edge.list)) %>%
      group_by(across(everything())) %>%
      summarise(n = n() / 3, .groups = "drop")
    colnames(adj.edge.list) = paste0("V", 1:4)
    networkData[setDT(as.data.frame(triangle.list)), on = paste0("V", 1:3), "h1" := i.V4]
    networkData[setDT(adj.edge.list), on = paste0("V", 1:3), "h2" := i.V4]
    networkData[is.na(networkData[["h1"]]), "h1"] = 0
    networkData[is.na(networkData[["h2"]]), "h2"] = 0
  }
  
  colnames(networkData)[1:3] = c("i", "j", "k")
  
  return(networkData)
}

### UGEE
ugee_triangle = function(formula,
                         data,
                         n,
                         id,
                         corstr = "independence",
                         zcor = NULL,
                         maxit = 200,
                         tol = 0.001,
                         true_param = NULL) {
  # Reorder network data by id
  network.data = setDT(as.data.frame(data))
  
  formula.h = unlist(strsplit(gsub(" ", "", as.character(formula[2]), fixed = T), "(?<=\\w)(?=\\W)|(?<=\\W)(?=\\w)", perl = T))
  formula.m = formula[-2]
  if (any(formula.h == "*") & length(formula.h) == 3) {
    if (any(colnames(network.data) == formula.h[1])) {
      formula.hh = colnames(network.data)[grepl(formula.h[1], colnames(network.data), fixed = T)]
      formula.hw = colnames(network.data)[grepl(formula.h[3], colnames(network.data), fixed = T)]
    } else {
      formula.hh = colnames(network.data)[grepl(formula.h[3], colnames(network.data), fixed = T)]
      formula.hw = colnames(network.data)[grepl(formula.h[1], colnames(network.data), fixed = T)]
    }
    mr1 = as.matrix(network.data[, ..formula.hh])[, 1]
    mr2 = as.matrix(network.data[, ..formula.hh])[, 2]
    if (sum(mr1) > sum(mr2)) {
      mr1 = as.matrix(network.data[, ..formula.hh])[, 2]
      mr2 = as.matrix(network.data[, ..formula.hh])[, 1]
    }
    mw = as.matrix(network.data[, ..formula.hw])
  } else if (length(formula.h) == 1) {
    formula.hh = colnames(network.data)[grepl(formula.h, colnames(network.data), fixed = T)]
    mr1 = as.matrix(network.data[, ..formula.hh])[, 1]
    mr2 = as.matrix(network.data[, ..formula.hh])[, 2]
    if (sum(mr1) > sum(mr2)) {
      mr1 = as.matrix(network.data[, ..formula.hh])[, 2]
      mr2 = as.matrix(network.data[, ..formula.hh])[, 1]
    }
    mw = matrix(nrow = 0, ncol = 0)
  } else {
    stop("Please enter the formula with the correct format!")
  }
  
  mm = model.matrix(formula.m, network.data)
  comb.index = as.matrix(network.data[, ..id])
  n.time = ncol(mw) + 1
  R = diag(3*n.time-1)
  if (n.time != 1) {
    beta10 = matrix(head(true_param, n.time))
    beta20 = matrix(true_param[(n.time+1):(2*n.time)])
    pi0 = tail(true_param, n.time)
  } else {
    beta10 = matrix(true_param[1])
    beta20 = matrix(true_param[2])
    pi0 = numeric()
  }
  
  # Estimate beta using Gauss-Newton method
  restart = 1
  while (restart == 1) {
    # Initialization
    if (!is.null(true_param)) {
      beta1 = beta10
      beta2 = beta20
      pi = pi0[-n.time]
    } else {
      beta1 = matrix(rnorm(n.time))
      beta2 = matrix(rnorm(n.time))
      pi = runif(n.time-1)
    }
    
    if (corstr != "independence") {
      alpha.matrix = cor_alpha(mr1, mr2, mw, mm, beta1, beta2, pi, n.time)
      if (corstr == "exchangeable") {
        alpha = sum(alpha.matrix) / (3*n.time-1) / (3*n.time-2)
        R[lower.tri(R)] = R[upper.tri(R)] = alpha
      } else if (corstr == "ar1") {
        alpha = sum(alpha.matrix) / (3*n.time-1) / (3*n.time-2)
        R[upper.tri(R)] = alpha^abs(row(R)[upper.tri(R)] - col(R)[upper.tri(R)])
        R[lower.tri(R)] = t(R)[lower.tri(R)]
      } else if (corstr == "unstructured") {
        R = R + alpha.matrix
      } else if (corstr == "fixed") {
        if (is.null(zcor)) stop("zcor must be provided for fixed correlation")
        R = zcor
      } else stop("correlation structure undefined!!")
    }
    
    # Main
    for (iter in 1:maxit) {
      asym.mean.result = Gauss_newton_triangle(mr1, mr2, mw, mm, beta1, beta2, pi, n.time, R)
      beta1.update = asym.mean.result$beta1
      beta2.update = asym.mean.result$beta2
      pi.update = asym.mean.result$pi
      
      R.update = R
      if (corstr == "exchangeable") {
        alpha.matrix = cor_alpha(mr1, mr2, mw, mm, beta1.update, beta2.update, pi.update, n.time)
        alpha = sum(alpha.matrix) / (3*n.time-1) / (3*n.time-2)
        R.update[lower.tri(R.update)] = R.update[upper.tri(R.update)] = alpha
      } else if (corstr == "ar1") {
        alpha.matrix = cor_alpha(mr1, mr2, mw, mm, beta1.update, beta2.update, pi.update, n.time)
        alpha = sum(alpha.matrix) / (3*n.time-1) / (3*n.time-2)
        R.update[upper.tri(R.update)] = alpha^abs(row(R.update)[upper.tri(R.update)] - col(R.update)[upper.tri(R.update)])
        R.update[lower.tri(R.update)] = t(R.update)[lower.tri(R.update)]
      } else if (corstr == "unstructured") {
        alpha.matrix = cor_alpha(mr1, mr2, mw, mm, beta1.update, beta2.update, pi.update, n.time)
        R.update = diag(3*n.time-1) + alpha.matrix
      }
      
      if (identical(beta1.update, beta1) & identical(beta2.update, beta2) & identical(pi.update, pi) == T) {
        restart = 1
        break
      } else if (max(norm(beta1.update - beta1, type = "I"), norm(beta2.update - beta2, type = "I"), norm(matrix(pi.update - pi), type = "I")) < tol) {
        restart = 0
        break
      } else if (iter < maxit) {
        restart = 0
        beta1 = beta1.update
        beta2 = beta2.update
        pi = pi.update
        R = R.update
      } else {
        restart = 0
        warning("Gauss-Newton does not converge!! Please increase number of iterations")
      }
    }
  }
  
  # Estimate asymptotic variance
  asym.var.result = asym_var_triangle(mr1, mr2, mw, mm, beta1.update, beta2.update, pi.update, n, n.time, comb.index, R.update)
  
  if (n.time != 1) {
    return(list(beta1 = beta1.update, beta2 = beta2.update, pi = pi.update, vtau = asym.var.result))
  } else {
    return(list(beta1 = beta1.update, beta2 = beta2.update, vbeta = asym.var.result))
  }
}

cor_alpha = function(mr1, mr2, mw, mm, beta1, beta2, pi, n.time) {
  if (n.time != 1) {
    h = cbind(floor(mw) * mr1, mr1, floor(mw) * mr2, mr2, mw)
    f = cbind(expit(mm %*% t(beta1)) * expit(mm %*% t(beta2)) * c(pi, 1)^3, expit(mm %*% t(beta2)) * c(pi, 1)^3, matrix(pi, nrow = nrow(mm), ncol = 2, byrow = T))
  } else {
    h = cbind(mr1, mr2)
    f = cbind(expit(mm %*% beta1) * expit(mm %*% beta2), expit(mm %*% beta2))
  }
  
  r = (h - f) / sqrt(f * (1 - f))
  r_clipped = pmax(pmin(r, 1), -1)
  alpha.matrix = alpha_estimation(r_clipped)
  
  return(alpha.matrix)
}

