# Helper functions for checkResidualVar

# Helper function to make S_gamma and its inverse -------------------------

makeSgamma_plots <- function(gamma, W, no_inv = F) {
  
  S_gamma <- diag(c(exp(W %*% gamma)))
  # S_gamma <- diag(exp(W %*% gamma))
  
  if(no_inv == F){
    to_return <- list(S_gamma = S_gamma)
  } else{
    
    S_gamma_inv <- diag(c(-exp(W %*% gamma)))
    # S_gamma_inv <- diag(-exp(W %*% gamma))
    to_return <- list(S_gamma = S_gamma, S_gamma_inv = S_gamma_inv) 
    
  }
  
  return(to_return)
}

# Helper function to make K at new predictor values ------------------------

makeK_plots <- function(r_val, Z1, Z2 = NULL, Znew_ind = 0, nugget = nugget,
                        get_inv = T){

  Z1r <- sweep(Z1, 2, sqrt(r_val), "*")
  if (is.null(Z2)) {
    Z2r <- Z1r
  } else {
    Z2r <- sweep(Z2, 2, sqrt(r_val), "*")
  }
  Kpart <- fields::rdist(Z1r, Z2r)^2
  K <- exp(-Kpart)
  
  if(get_inv == T){
    
    if(is.null(Z2) == T & Znew_ind == 0){
      K <- K + nugget*diag(nrow(K)) 
      chol_K <- chol(K)
      K_inv <- chol2inv(chol_K)
      K_log_det <- 2*sum(log(diag(chol_K)))
      
      to_return <- list(K = K, K_inv = K_inv, K_log_det = K_log_det)
    } else{
      to_return <- list(K = K)
    }
    
  } else{
    
    to_return <- list(K = K)
  }  
  
  return(to_return)
}