#' Modification of the function "SimData" from the BKMR package to simulate heteroscedastic response data
#'
#' This function simulates a dataset with specified parameters for a given family 
#' and distribution type. 
#'
#' @param n Integer. The number of observations (default is 100).
#' @param M Integer. The number of mixture components (default is 5).
#' @param beta.true Numeric. The true value of the regression coefficient (default is 2).
#' @param hfun Integer. The functional form of the h function data should be generated under (default is 3).
#' @param Zgen Character. Determines how the matrix of mixtures is generated. It can be 
#'        one of the following strings: "unif" (uniform distribution), "norm" (normal distribution), 
#'        "corr" (correlated variables), or "realistic" (realistic mixture generation) (default is "norm").
#' @param ind Numeric vector. The indices of covariates h is a function of (default is 1:2).
#' @param family Character. The type of distribution used for generating covariates, 
#'        with options "norm" (normal), and "binom" (binomial), (default is "norm").
#' @param var_metal Numeric. The metal concentration, which determines the variance 
#'        as a function of the concentration (default is 3).
#' @param var_mult Numeric. A multiplier on the metal concentration to compute 
#'        the variance (default is 0.2).
#'
#' @return A simulated dataset with specified properties, including outcome variables 
#'         and covariates.
#'
#' @examples
#' # Generate a synthetic dataset with default parameters
#' data <- SimHData(n = 200)
#'
#' # Generate a dataset with a realistic mixture matrix and simple h function
#' data <- SimHData(n = 200, Zgen = "realistic, hfun = 1)
#'
SimHData <- function (n = 100, M = 5, beta.true = 2, hfun = 3, 
                      Zgen = "norm", ind = 1:2, family = "gaussian", var_metal = 3,
                      var_mult = 0.2) 
{
  
  HFun1 <- function(z, ind = 1) 4*plogis(z[ind[1]], 0, 0.3)
  HFun2 <- function(z, ind = 1:2) 1/4*(z[ind[1]] + z[ind[2]] + 1/2*z[ind[1]]*z[ind[2]])
  HFun3 <- function(z, ind = 1:2) 4*plogis(1/4*(z[ind[1]] + z[ind[2]] + 1/2*z[ind[1]]*z[ind[2]]), 0, 0.3)
  
  
  stopifnot(n > 0, M > 0, family %in% c("gaussian", "binomial"))
  if (family == "binomial") {
    sigsq.true <- 1
  }
  if (hfun == 1) {
    HFun <- HFun1
  }
  else if (hfun == 2) {
    HFun <- HFun2
  }
  else if (hfun == 3) {
    HFun <- HFun3
  }
  else {
    stop("hfun must be an integer from 1 to 3")
  }
  if (Zgen == "unif") {
    Z <- matrix(runif(n * M, -2, 2), n, M)
  }
  else if (Zgen == "norm") {
    Z <- matrix(rnorm(n * M), n, M)
  }
  else if (Zgen == "corr") {
    if (M < 3) {
      stop("M must be an integer > 2 for Zgen = 'corr'")
    }
    Sigma <- diag(1, M, M)
    Sigma[1, 3] <- Sigma[3, 1] <- 0.95
    Sigma[2, 3] <- Sigma[3, 2] <- 0.3
    Sigma[1, 2] <- Sigma[2, 1] <- 0.1
    Z <- MASS::mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  }
  else if (Zgen == "realistic") {
    VarRealistic <- structure(c(0.72, 0.65, 0.45, 0.48, 0.08, 
                                0.14, 0.16, 0.42, 0.2, 0.11, 0.35, 0.1, 0.11, 0.65, 
                                0.78, 0.48, 0.55, 0.06, 0.09, 0.17, 0.2, 0.16, 0.11, 
                                0.32, 0.12, 0.12, 0.45, 0.48, 0.56, 0.43, 0.11, 0.15, 
                                0.23, 0.25, 0.28, 0.16, 0.31, 0.15, 0.14, 0.48, 0.55, 
                                0.43, 0.71, 0.2, 0.23, 0.32, 0.22, 0.29, 0.14, 0.3, 
                                0.22, 0.18, 0.08, 0.06, 0.11, 0.2, 0.95, 0.7, 0.45, 
                                0.22, 0.29, 0.16, 0.24, 0.2, 0.13, 0.14, 0.09, 0.15, 
                                0.23, 0.7, 0.8, 0.36, 0.3, 0.35, 0.13, 0.23, 0.17, 
                                0.1, 0.16, 0.17, 0.23, 0.32, 0.45, 0.36, 0.83, 0.24, 
                                0.37, 0.2, 0.36, 0.34, 0.25, 0.42, 0.2, 0.25, 0.22, 
                                0.22, 0.3, 0.24, 1.03, 0.41, 0.13, 0.39, 0.1, 0.1, 
                                0.2, 0.16, 0.28, 0.29, 0.29, 0.35, 0.37, 0.41, 0.65, 
                                0.18, 0.3, 0.18, 0.16, 0.11, 0.11, 0.16, 0.14, 0.16, 
                                0.13, 0.2, 0.13, 0.18, 0.6, 0.18, 0.13, 0.08, 0.35, 
                                0.32, 0.31, 0.3, 0.24, 0.23, 0.36, 0.39, 0.3, 0.18, 
                                0.79, 0.42, 0.12, 0.1, 0.12, 0.15, 0.22, 0.2, 0.17, 
                                0.34, 0.1, 0.18, 0.13, 0.42, 1.27, 0.1, 0.11, 0.12, 
                                0.14, 0.18, 0.13, 0.1, 0.25, 0.1, 0.16, 0.08, 0.12, 
                                0.1, 0.67), .Dim = c(13L, 13L))
    if (M > ncol(VarRealistic)) {
      stop("Currently can only generate exposure data based on a realistic correlation structure with M = 13 or fewer. Please set M = 13 or use Zgen = c('unif','norm'")
    }
    else if (M <= 13) {
      Sigma <- VarRealistic[1:M, 1:M]
    }
    Z <- MASS::mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  }
  colnames(Z) <- paste0("z", 1:M)
  X <- cbind(3 * cos(Z[, 1]) + 2 * rnorm(n))
  
  sigma_2 <- exp(1 + var_mult*Z[,var_metal]) 
  eps <- rnorm(n, sd = sqrt(sigma_2))
  
  h <- apply(Z, 1, HFun)
  mu <- X * beta.true + h
  y <- drop(mu + eps)
  if (family == "binomial") {
    ystar <- y
    y <- ifelse(ystar > 0, 1, 0)
  }
  dat <- list(n = n, M = M, sigsq.true = sigma_2, beta.true = beta.true, 
              Z = Z, h = h, X = X, y = y, hfun = hfun, HFun = HFun, 
              family = family)
  if (family == "binomial") {
    dat$ystar <- ystar
  }
  dat
}

