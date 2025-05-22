
makePredictionInvervalPlot <- function(beta, gamma, r, tau, y, X, Z, W,
                                                        sel = NULL, nugget = 0) {
  
  Z_new <- Z
  X_new <- X
  W_new <- W
  
  # This is to store the element-wise summations across the for loop 
  sum_var_covar_mat <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  
  n_samps <- length(tau)
  
  if(is.null(sel)){
    sel_index <- (1:n_samps)[seq(1, n_samps, by = 10)]
    warning("sel = NULL. Defaulting to every 10th posterior sample for exact computation.")
  } else {
    sel_index <- sel
  }
  
  n_thin_samps <- length(sel_index)
  
  post_comps_store <- list(postmean = vector("list", n_thin_samps))
  
  for (s in sel_index) {
    K_znew_z <- makeK_plots(r[s, ], Z_new, Z, nugget = nugget)$K
    K_z_znew <- makeK_plots(r[s, ], Z, Z_new, nugget = nugget)$K
    K_z <- makeK_plots(r[s, ], Z, nugget = nugget, Znew_ind = 1)$K # Same as K_z_new when the same Z is inputted
    
    S_gamma <- makeSgamma_plots(gamma[s,], W, no_inv = T)$S_gamma # Same as S_gamma_new when the same W is inputted
    
    V <- tau[s]*K_z + S_gamma
    V_inv <- chol2inv(chol(V))
    
    postvar <- as.matrix(tau[s] * K_z + S_gamma - tau[s]^2 * K_znew_z %*% V_inv %*% K_z_znew)
    postmean <- as.vector(X_new %*% beta[s,] + tau[s]*K_znew_z %*% V_inv %*% (y - X %*% beta[s,]))
    
    sum_var_covar_mat <- sum_var_covar_mat + postvar
    
    post_comps_store$postmean[[s]] <- postmean
    
  }
  
  postmean <- post_comps_store$postmean[sel_index]
  
  postmean_mat <- t(do.call("cbind", postmean))
  m <- colMeans(postmean_mat)
  ve <- var(postmean_mat)
  ev <- sum_var_covar_mat / n_thin_samps
  v <- ve + ev
  ret <- data.frame(y, prediction = m, sd = sqrt(diag(v)))
  return(ret)
}



makeOverallRiskPlot <- function(beta, gamma, r, tau, y, X, Z, W, q.fixed = 0.5, q.levels = seq(0.1, 0.9, by = 0.1),
                                sel = NULL, method = "exact", nugget = 0) {
  
  cc <- c(1, -1)
  if (!{method %in% c("approx", "exact")}) {
    stop(stop("method must be one of c('approx', 'exact')"))
  }
  if (length(q.fixed) > 1) {
    stop("'q.fixed' must be one-dimensional.")
  }
  
  q_vals_fixed <- apply(Z, 2, function(x) quantile(x, probs = q.fixed))
  
  q_vals_levels <- apply(Z, 2, function(x) quantile(x, probs = q.levels))
  
  mu_samp <- matrix(NA, nrow = 2, ncol = length(tau))
  var_samp <- array(NA, dim = c(2, 2, length(tau)))
  
  df_list <- vector("list", length = length(q.levels))
  
  n_samps <- length(tau)
  
  if(is.null(sel)){
    sel_index <- (1:n_samps)[seq(1, n_samps, by = 10)]
    warning("sel = NULL. Defaulting to every 10th posterior sample for exact computation.")
  } else {
    sel_index <- sel
  }
  
  n_thin_samps <- length(sel_index)
  
  for (i in 1:nrow(q_vals_levels)) {
    Z_new <- rbind(varying_quantile = q_vals_levels[i, ], q_vals_fixed)
    if (method == "exact") {
      for (s in sel_index) {
        K_znew_z <- makeK_plots(r[s, ], Z_new, Z, nugget = nugget)$K
        K_z_znew <- makeK_plots(r[s, ], Z, Z_new, nugget = nugget)$K
        K_z <- makeK_plots(r[s, ], Z, nugget = nugget, Znew_ind = 1)$K
        K_znew <- makeK_plots(r[s, ], Z_new, Znew_ind = 1, nugget = nugget)$K
        
        V <- tau[s]*K_z + makeSgamma_plots(gamma[s,], W, no_inv = T)$S_gamma
        V_inv <- chol2inv(chol(V))
        
        var_samp[, , s] <- tau[s] * K_znew - tau[s]^2 * K_znew_z %*% V_inv %*% K_z_znew
        mu_samp[, s] <- tau[s]*K_znew_z %*% V_inv %*% (y - X %*% beta[s,])
        
      }
      
      est <- apply(mu_samp, 1, mean, na.rm = T)
      diff <- cc %*% est
      var_mean <- var(t(mu_samp), na.rm = T)
      
      mean_var <- apply(var_samp, c(1, 2), mean, na.rm = T)
      
      
      var <- var_mean + mean_var
      diff_var <- cc %*% var %*% cc
      diff_sd <- sqrt(diff_var) }
    
    message(paste0("Computing ", as.name(q.levels[i]), "% : ", i, " of ", length(q.levels)))
    
    df_list[[i]] <- data.frame(
      q.level = q.levels[i],
      est = diff, sd = diff_sd,
      var = diff_var
    )
  }
  
  return(do.call(rbind, df_list))
}

