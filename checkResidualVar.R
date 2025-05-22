#' Check Residual Variance for BKMR Model
#'
#' This function checks the residual variance from a BKMR model and provides residuals plots
#' for visual diagnostics based on specified model settings and residual calculation methods.
#'
#' @import np
#'
#' @param fit Output from NIMBLE::runMCMC
#' @param y A numeric vector representing the outcome variable.
#' @param X A matrix of confounders, with column names for each confounder.
#' @param Z A matrix of mixture components, with column names for each mixture component.
#' @param method A character string specifying the method for fitting the model, either 'exact' or 'approx' as in the BKMR package (default is "exact").
#' @param Z_component A numeric vector of indices specifying which columns of the mixture matrix to plot against the residuals (default is NULL).
#' @param X_component A numeric vector of indices specifying which columns of the covariate matrix to plot against the residuals on the x-axis (default is NULL).
#'        Set this to NULL if using the fitted values.
#' @param return_residuals Logical. If TRUE, a vector of residual values will be returned, otherwise FALSE (default is TRUE).
#' @param residuals A character string indicating which method to use for calculating residuals: 'bkmr' for BKMR residuals, 'kmr' for KMR residuals, and 'linreg' for linear regression residuals (default is "bkmr").
#' @param int Logical. If TRUE, a linear regression model with interactions is fit. If FALSE, a main effects model is fit (default is TRUE).
#' @param sel Indices to thin the posterior samples for faster exact computation (default is to keep every 10th sample).
#' @param categorical_X A list of names and indices of categorical variables that have been included as dummies in the X matrix (default is NULL).
#' @param ... Additional arguments passed to the KMR function.
#'
#' @return A plot based on the specified plot type and, optionally, a vector of residuals.
#'
#' @examples
#' # Example usage of the function
#' fit <- BKMR::bkmr() # Example BKMR model fit
#' residual_plot <- checkResidualVar(fit = fit, y = outcome, X = confounders, Z = mixtures, method = 'exact', 
#'                                   residuals = 'bkmr', return_residuals = TRUE)
#'
checkResidualVar <- function(fit = NULL, y, X, Z, method = "exact",
                             Z_component = NULL, X_component = NULL,
                             return_residuals = T, residuals = "bkmr", 
                             int = T, sel = NULL, categorical_X = NULL, ...){
  
  # Input checks
  
  if (is.null(fit) & residuals == "bkmr") {
    stop("fit object must be supplied if residuals = `bkmr`")
  }
  
  stopifnot("`y` must be a vector" = is.vector(y))
  
  stopifnot("`X` must be a matrix" = is.matrix(X))
  
  stopifnot("`Z` must be a matrix" = is.matrix(Z))
  
  if (!(method %in% c("exact", "approx"))) {
    warning("`method` not recognized, defaulting to `exact`.")
    method <- "exact"
  }
  
  if (is.null(method) & residuals == "bkmr") {
    warning("`method` must be supplied with residuals = `bkmr`. Defaulting to `exact`")
    method <- "exact"
  }
  
  stopifnot("`y`, `X`, and `Z` must have the same length." = length(y) == nrow(X) & nrow(X) == nrow(Z))
  
  if (!all(Z_component %in% seq_len(ncol(Z)))) {
    warning("Indices supplied in `Z_component` not in ncol(Z), defaulting to NULL")
    Z_component <- NULL
  }
  
  if (!all(X_component %in% seq_len(ncol(X)))) {
    warning("Indices supplied in `X_component` not in ncol(X), defaulting to NULL")
    X_component <- NULL
  }
  
  if (!(residuals %in% c("bkmr", "kmr", "linreg"))) {
    warning("Argument supplied to `residuals` is invalid, defaulting to `bkmr`")
  }
  
  if (!is.logical(return_residuals)) {
    warning("`return_residuals` must be logical, defaulting to TRUE")
    return_residuals <- TRUE
  }
  
  stopifnot("`X` must have column names." = !is.null(colnames(X)))
  
  stopifnot("`Z` must have column names." = !is.null(colnames(Z)))
  
  approx.bkmr <- function(){
    if(length(fit) == 2 & "WAIC" %in% names(fit)){
      # If WAIC was being monitored 
      ests <- fit$samples
    } else {
      # If WAIC was not being monitored 
      ests <- fit
    }
    
    # Make sure the sampler is monitoring sigma2
    if(all(isFALSE(grepl("sigma2", colnames(ests))))){
      warning("Posterior samples for `sigma2` are missing. Ensure NIMBLE monitors `sigma2`")
    }
    
    sigsq_eps <- mean(ests[, grepl("sigma2", colnames(ests))][,1])
    
    r <- colMeans(ests[, grepl("r", colnames(ests))])
    beta <- colMeans(ests[, grepl("beta", colnames(ests))])
    tau <- mean(ests[, grepl("tau", colnames(ests))])
    gamma <- mean(ests[, grepl("gamma", colnames(ests))])
    lambda <- mean(ests[, grepl("tau", colnames(ests))] / 
                     ests[, grepl("sigma2", colnames(ests))][,1])
    
    # The outcome vector has to be supplied as an argument 
    ycont <- y
    K <- makeK_plots(r, Z1 = Z, Z2 = NULL, Znew_ind = 0, nugget = 0, get_inv = F)$K
    
    V <- diag(1, nrow(Z), nrow(Z)) + lambda * K
    
    # V <- sigsq_eps * diag(1, nrow(Z), nrow(Z)) + tau * K (defined differently by BKMR authors)
    cholV <- chol(V)
    Vinv <- chol2inv(cholV)
    lamKVinv <- lambda*K%*%Vinv
    postvar <- lambda*sigsq_eps*(K - lamKVinv%*%K)
    postmean <- lamKVinv %*% (ycont - X%*%beta)
    ret <- list(postmean = drop(postmean), postvar = postvar)
    ret
    
  }
  exact.bkmr <- function(){
    
    if(length(fit) == 2 & "WAIC" %in% names(fit)){
      # If WAIC was being monitored 
      ests <- fit$samples
    } else {
      # If WAIC was not being monitored 
      ests <- fit
    }
    
    # Make sure the sampler is monitoring sigma2
    if(all(isFALSE(grepl("sigma2", colnames(ests))))){
      warning("Posterior samples for `sigma2` are missing. Ensure NIMBLE monitors `sigma2`")
    }
    
    n_samps <- nrow(fit$samples)
    
    if(is.null(sel)){
      sel_index <- (1:n_samps)[seq(1, n_samps, by = 10)]
      warning("sel = NULL. Defaulting to every 10th posterior sample for exact computation.")
    } else {
      sel_index <- sel
    }
    
    n_thin_samps <- length(sel_index)
    
    post_comps_store <- list(postmean = vector("list", n_thin_samps),
                             postvar = vector("list", n_thin_samps))
    
    sigsq_eps <- ests[, grepl("sigma2", colnames(ests))][,1]
    r <- ests[, grepl("r", colnames(ests))]
    beta <- ests[, grepl("beta", colnames(ests))]
    tau <- ests[, grepl("tau", colnames(ests))]
    gamma <- ests[, grepl("gamma", colnames(ests))]
    lambda <- ests[, grepl("tau", colnames(ests))] / sigsq_eps
    ycont <- y
    
    for(i in sel_index){
      K <- makeK_plots(r_val = r[i, ], Z1 = Z, Z2 = NULL, Znew_ind = 0, nugget = 0,
                       get_inv = F)$K
      V <- diag(1, nrow(Z), nrow(Z)) + lambda[i] * K
      cholV <- chol(V)
      Vinv <- chol2inv(cholV)
      lamKVinv <- lambda[i]*K%*%Vinv
      postvar <- lambda[i]*sigsq_eps[i]*(K - lamKVinv%*%K)
      postmean <- lamKVinv %*% (ycont - X%*%beta[i, ])
      post_comps_store$postmean[[i]] <- postmean
      post_comps_store$postvar[[i]] <- postvar
    }
    
    postmean <- post_comps_store$postmean[sel_index]
    postvar <- post_comps_store$postvar[sel_index]
    
    postmean_mat <- t(do.call("cbind", postmean))
    m <- colMeans(postmean_mat)
    ve <- var(postmean_mat)
    # ev <- apply(postvar_arr, c(1, 2), mean)
    # This is a more efficient way to compute the element-wise mean 
    ev <- Reduce("+", postvar) / length(postvar)
    v <- ve + ev
    ret <- list(postmean = m, postvar = v)
    
    ret
  }
  
  bkmr.residuals <- function(){
    
    if (method == "exact") {
      h_mean <- exact.bkmr()$postmean
    } else {
      h_mean <- approx.bkmr()$postmean
    }
    
    
    ests <- fit$samples
    beta_mean <- colMeans(ests[, grepl("beta", colnames(ests))])
    Xbeta_mean <- X %*% beta_mean 
    
    res <- y - (h_mean + Xbeta_mean)
    fit <- h_mean + Xbeta_mean
    
    df <- data.frame(y, X_aug, Z, res, fit)
    
    bkmr.plots.scatterplot(df)
  }
  kmr.residuals <- function(){
    
    kmr_bandwidth <- npregbw(xdat = data.frame(Z, X), ydat = y, ...)
    
    kmr_output <- npreg(kmr_bandwidth, txdat = data.frame(Z, X), tydat = y, residuals = TRUE)
    
    res <- kmr_output$resid
    
    fit <- fitted(kmr_output)
    
    df <- data.frame(y, X_aug, Z, res, fit)
    
    kmr.plots.scatterplot(df)
  }
  lm.residuals <- function(){
    design_matrix <- model.matrix(~ .^2, as.data.frame(Z))
    predictors_Z <- design_matrix[, -1]
    
    if (int == T){
      linreg <- lm(y ~ predictors_Z + X)
    } else if (int == F) {
      linreg <- lm(y ~ X + Z)
    }
    
    res <- linreg$residuals
    fit <- linreg$fitted.values
    
    df <- data.frame(y, X_aug, Z, res, fit)
    
    lm.plots.scatterplot(df)
  }
  
  bkmr.plots.scatterplot <- function(df){
    fitted_plot <- local({
      fit_plot <- ggplot(df, aes(x = fit, y = res)) +
        geom_point() +
        xlab("Posterior mean of y") +
        ylab("Bayesian residual")
      
      print(fit_plot)
    })
    
    if (num_mixture_plots > 0) {
      for (i in 1:num_mixture_plots) {
        mixture_plots[[i]] <- local({
          i <- i
          
          mix_plot <- ggplot(df, aes(x = get(mixture_names[i]), y = res)) +
            geom_point() +
            xlab(mixture_names[i]) +
            ylab("Bayesian residual")
          
          print(mix_plot)
        })
      }
    }
    
    if (num_covariate_plots > 0) {
      for (i in 1:num_covariate_plots) {
        covariate_plots[[i]] <- local({
          i <- i
          
          if(covariate_names[i] %in% categorical_covariate_names){
            
            covariate_plot <- ggplot(df, aes(x = factor(get(covariate_names[i])), y = res)) +
              geom_boxplot() +
              xlab(covariate_names[i]) +
              ylab("Bayesian residual")
            
            
            print(covariate_plot)
            
          } else {
            
            covariate_plot <- ggplot(df, aes(x = get(covariate_names[i]), y = res)) +
              geom_point() +
              xlab(covariate_names[i]) +
              ylab("Bayesian residual")
            
            print(covariate_plot)
          }
          
          
        })
      }
    }
    
    if (num_mixture_plots == 0 & num_covariate_plots == 0) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, bayes_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot)
      }
    }
    
    if (num_mixture_plots > 0 & num_covariate_plots > 0) {
      if (return_residuals == T) {
        to_return <- list(
          fitted_plot = fitted_plot, mixture_plots = mixture_plots, covariate_plots = covariate_plots,
          bayes_res = df$res
        )
      } else {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots, covariate_plots = covariate_plots)
      }
    } else if (num_mixture_plots > 0 & (num_covariate_plots == 0 | is.null(num_covariate_plots))) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots, bayes_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots)
      }
    } else if ((num_mixture_plots == 0 | is.null(num_mixture_plots)) & num_covariate_plots > 0) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, covariate_plots = covariate_plots, bayes_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot, covariate_plots = covariate_plots)
      }
    }
    
    return(to_return)
  }
  kmr.plots.scatterplot <- function(df){
    
    fitted_plot <- local({
      fit_plot <- ggplot(df, aes(x = fit, y = res)) +
        geom_point() +
        xlab("Fitted") +
        ylab("Residual")
      
      print(fit_plot)
    })
    
    if (num_mixture_plots > 0) {
      for (i in 1:num_mixture_plots) {
        mixture_plots[[i]] <- local({
          i <- i
          
          mix_plot <- ggplot(df, aes(x = get(mixture_names[i]), y = res)) +
            geom_point() +
            xlab(mixture_names[i]) +
            ylab("Residual")
          
          print(mix_plot)
        })
      }
    }
    
    if (num_covariate_plots > 0) {
      for (i in 1:num_covariate_plots) {
        covariate_plots[[i]] <- local({
          i <- i
          
          if(covariate_names[i] %in% categorical_covariate_names){
            
            covariate_plot <- ggplot(df, aes(x = factor(get(covariate_names[i])), y = res)) +
              geom_boxplot() +
              xlab(covariate_names[i]) +
              ylab("Residual")
            
            
            print(covariate_plot)
            
          } else {
            
            covariate_plot <- ggplot(df, aes(x = get(covariate_names[i]), y = res)) +
              geom_point() +
              xlab(covariate_names[i]) +
              ylab("Residual")
            
            print(covariate_plot)
          }
          
          
        })
      }
    }
    if (num_mixture_plots == 0 & num_covariate_plots == 0) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, kmr_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot)
      }
    }
    
    if (num_mixture_plots == 0 & num_covariate_plots == 0) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, kmr_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot)
      }
    }
    
    if (num_mixture_plots > 0 & num_covariate_plots > 0) {
      if (return_residuals == T) {
        to_return <- list(
          fitted_plot = fitted_plot, mixture_plots = mixture_plots, covariate_plots = covariate_plots,
          kmr_res = df$res
        )
      } else {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots, covariate_plots = covariate_plots)
      }
    } else if (num_mixture_plots > 0 & (num_covariate_plots == 0 | is.null(num_covariate_plots))) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots, kmr_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots)
      }
    } else if ((num_mixture_plots == 0 | is.null(num_mixture_plots)) & num_covariate_plots > 0) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, covariate_plots = covariate_plots, kmr_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot, covariate_plots = covariate_plots)
      }
    }
    return(to_return)
  }
  lm.plots.scatterplot <- function(df){
    fitted_plot <- local({
      fit_plot <- ggplot(df, aes(x = fit, y = res)) +
        geom_point() +
        xlab("Fitted") +
        ylab("Residual")
      
      print(fit_plot)
    })
    
    if (num_mixture_plots > 0) {
      for (i in 1:num_mixture_plots) {
        mixture_plots[[i]] <- local({
          i <- i
          
          mix_plot <- ggplot(df, aes(x = get(mixture_names[i]), y = res)) +
            geom_point() +
            xlab(mixture_names[i]) +
            ylab("Residual")
          
          print(mix_plot)
        })
      }
    }
    
    if (num_covariate_plots > 0) {
      for (i in 1:num_covariate_plots) {
        covariate_plots[[i]] <- local({
          i <- i
          
          if(covariate_names[i] %in% categorical_covariate_names){
            
            covariate_plot <- ggplot(df, aes(x = factor(get(covariate_names[i])), y = res)) +
              geom_boxplot() +
              xlab(covariate_names[i]) +
              ylab("Residual")
            
            
            print(covariate_plot)
            
          } else {
            
            covariate_plot <- ggplot(df, aes(x = get(covariate_names[i]), y = res)) +
              geom_point() +
              xlab(covariate_names[i]) +
              ylab("Residual")
            
            print(covariate_plot)
          }
          
          
        })
      }
    }
    if (num_mixture_plots == 0 & num_covariate_plots == 0) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, linreg_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot)
      }
    }
    
    if (num_mixture_plots > 0 & num_covariate_plots > 0) {
      if (return_residuals == T) {
        to_return <- list(
          fitted_plot = fitted_plot, mixture_plots = mixture_plots, covariate_plots = covariate_plots,
          linreg_res = df$res
        )
      } else {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots, covariate_plots = covariate_plots)
      }
    } else if (num_mixture_plots > 0 & (num_covariate_plots == 0 | is.null(num_covariate_plots))) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots, linreg_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot, mixture_plots = mixture_plots)
      }
    } else if ((num_mixture_plots == 0 | is.null(num_mixture_plots)) & num_covariate_plots > 0) {
      if (return_residuals == T) {
        to_return <- list(fitted_plot = fitted_plot, covariate_plots = covariate_plots, linreg_res = df$res)
      } else {
        to_return <- list(fitted_plot = fitted_plot, covariate_plots = covariate_plots)
      }
    }
    
    return(to_return)
  }
  
  num_mixture_plots <- length(Z_component)
  mixture_names <- colnames(Z[, Z_component, drop = F])
  
  num_rows <- nrow(X)  
  num_cols <- length(categorical_X)  
  cat_df <- data.frame(matrix(NA_character_, nrow = num_rows, ncol = num_cols))
  for(i in 1:length(categorical_X)){
    # Extract indices of categorical columns
    indices <- categorical_X[[i]]  
    categorical_covariates <- X[, indices, drop = FALSE]  
    
    # Re-construct the reference level column from the dummy-encoded variables
    ref_dummy <- apply(categorical_covariates, 1, function(x) all(x == 0))  
    joint_cat <- cbind(ref = ifelse(ref_dummy, 1, 0), categorical_covariates)
    
    # Convert each row to the corresponding category name
    cat_vec <- apply(joint_cat, 1, function(x) colnames(joint_cat)[which(x == 1)]) |> as.factor()
    
    # Assign to the output dataframe
    cat_df[, i] <- cat_vec
    colnames(cat_df)[i] <- names(categorical_X)[i]
  }
  
  # Indices of the X matrix supplied as being categorical 
  cat_indices <- unname(unlist(categorical_X))
  # Indices supplied in X_component but not in categorical_X are assumed to be cont.
  cont_indices <- X_component[!(X_component %in% cat_indices)]
  
  cont_df <- data.frame(X[, cont_indices, drop = F])
  
  # Final X object based on X_components and categorical_X for plotting purposes.
  X_aug <- cbind(cat_df, cont_df)
  
  num_covariate_plots <- ncol(X_aug)
  covariate_names <- colnames(X_aug)
  
  categorical_covariate_names <- names(categorical_X)
  
  mixture_plots <- list()
  covariate_plots <- list()
  
  if(residuals == "bkmr"){
    bkmr.residuals()
  } else if(residuals == "kmr"){
    kmr.residuals()
  } else if(residuals == "linreg"){
    lm.residuals()
  }
  
}
