rm(list = ls())

library(MASS)
library(Matrix)
library(extraDistr)
library(expm)
library(SimCorMultRes)

#set.seed(10)

################################################################
#########  normal and logistic simulation
################################################################
#' @param n_sub number of subject
#' if using the sum to 0 constraint. Ignored if constrain = FALSE.
#' @param min_obs minimum number of observations per subject
#' @param max_obs maximum number of observations per subject
#' @param n_pred number of predictors (p)
#' @param sigma_X marginal sd for the predictors
#' @param sigma_Y marginal sd for the outcome  
#' @param corr_X correlation structure for predictors
#' should be a p*p matrix, if not specified, will use identity
#' @param beta coefficients for the predictors
#' @param alpha_obs correlation parameter for observations of the 
#' same subject, only support exchangeable (or independence)
#' @param link 'identity' for normal and 'logit' for binary

sim_norm_logit <- function (n_sub = 50, min_obs = 1, max_obs = 6,
                         n_pred = 40, sigma_X = 1, sigma_Y = 2, corr_X = NA,
                         alpha_obs = 0.5, link = 'identity',
                         beta = c(rep(0, times = 10), rep(-2, times = 10),
                                  rep(0, times = 10), rep(2, times = 10))) {
  
  # create placeholders
  beta <- matrix(beta, nrow = n_pred, ncol = 1)
  id_vec <- vector()
  measure_idx_vec <-vector()
  err_vec <- vector()
  X <- matrix(, nrow = 0, ncol = n_pred)
  Y_binom <- vector()
  
  #if not specify correlation for X, use independence
  if (is.na(corr_X) ) {
    corr_X <- diag(rep(1, times = length(beta)))
  }
  
  # generate data by subject
  
  for (id in 1: n_sub) {
    
    # number of observations for subject i
    ni <- rdunif(1, min_obs, max_obs)
    
    # id index
    id_vec <- c(id_vec,
                rep(id, each = ni))
    
    # measurement index
    measure_idx_vec <- c(measure_idx_vec,
                         c(1: ni))
    
    # generate correlated predictors
    X_i <- mvrnorm(n = ni,
                   mu = rep(0, times = n_pred),
                   Sigma = corr_X * sigma_X ^ 2)
    
    X <- rbind(X, X_i)
    
    # correlation between observations
    R <- matrix(alpha_obs, nrow = ni, ncol = ni)
    diag(R) <- 1
    
    # error for normal data 
    if (link == 'identity') {
      err <- mvrnorm(n = 1,
                     mu = rep(0, times = ni),
                     Sigma = R * sigma_Y ^ 2)
      err <- as.vector(t(err))
      err_vec <- c(err_vec, err)
    }
    
    # random outcome for correlated binary data
    
    if (link == 'logit') {
      ## single observation
      if (ni == 1) {
        pi <- exp(X_i %*% beta)/(1+exp(X_i %*% beta))
        Y_i <- rbinom(n = 1, size = 1, prob = pi)
        Y_binom <- cbind(Y_binom, Y_i)
        
      } else {
        ## multiple observations
        sim_bin <- rbin(clsize = ni, intercepts = 0, betas = as.vector(beta / 4),
                        xformula = ~X_i,
                        cor.matrix = R, link = "logit")
        Y_i <- sim_bin$Ysim
        Y_binom <- cbind(Y_binom, Y_i)
        
      }
    }
  }
  
  SNR <- sqrt(t(X %*% beta) %*% X %*% beta / (length(id_vec) * sigma_Y ^ 2))
  
  if (link == 'identity') {
    # generate outcome Y
    Y <- X %*% beta + err_vec
    # calculate signal - to - noise ratio
  }
  
  if (link == 'logit') {
    # generate outcome Y
    Y <- as.vector(Y_binom)
  }
  
  # prepare aggreagated dataset
  long_dat <- data.frame(
    id = id_vec,
    Y = Y,
    measure_idx = measure_idx_vec)
  
  long_dat <- cbind(long_dat, X)
  colnames(long_dat) <- c('id', 'Y', 'measure_idx',
                          paste("x", 1: n_pred, sep = ""))
  rownames(long_dat) <- NULL
  
  if (link == 'identity') {
    long_dat$error <- err_vec
  }  
  
  if (link == 'logit') {
    long_dat$margin_p <- exp(X %*% beta)/(1+exp(X %*% beta))
  }
  
  # summarize
  corr_Y <- 'exchangeable'
  if (alpha_obs == 0) {
    corr_Y <- 'independence'
  }
  
  cat("simulated", n_sub, "subjects, ")
  cat("in total there are", length(Y), "observations", "\n")
  cat("used", corr_Y, "correlation between observations ")
  if (alpha_obs != 0) {
    cat("with alpha =",
        alpha_obs)
  }
  cat("\n")
  if(link=='identity'){cat("signal-to-noise ratio is ", SNR)}
  
  return (list(
    X = X,
    Y = Y,
    id = id_vec,
    SNR = SNR,
    whole_dat = long_dat,
    beta = beta
  ))
  
}



######## Example
tmp_logit <- sim_norm_logit(link='logit',alpha_obs = 0.5)
ldat<-tmp_logit$whole_dat

if (FALSE) {
  source('../pc_gee.R')
  
  res_tmp <- fit.gee(y = tmp_logit$Y, X = tmp_logit$X, 
                     id = tmp_logit$id,
                     constrain = TRUE, rho = 100,link='logit',
                     gamma = 1, lambda = 0.25, corstr = 'independence')
  sum(res_tmp$beta)
}

################################################################
#########  correlation structure for X
################################################################
#' @param n_pred number of predictors
#' @param alpha correlation parameter in exchangeable and AR-1
#' @param corstr correlation structure
#' allow 'independence', 'AR-1' and 'block diagonal '
#' @param blocksize a vector of block sizes for block diagonal  

build_corr_X <- function (n_pred = 40, alpha = 0,
                           corstr = 'independence',
                           blocksize = NA) {
  
  ## block diagonal
  if (corstr == "block diagonal" && is.na(blocksize))
    stop("Need to provide block size for block diagnanol correlation")
  
  if (corstr == "block diagonal") {
    block_list <- list()
    for (j in 1: length(blocksize)) {
      block = blocksize[j]
      if (block == 1) {
        S_i = 1
      } else {
        S_i <- toeplitz(c(1, rep(alpha, times = block - 1)))
      }
      block_list[[j]] <- S_i
    }
    return (as.matrix(bdiag(block_list)))
  }
  
  ## AR-1
  if (corstr == 'AR-1') {
    ar1_cor <- function (n, rho) {
      exponent <- abs(matrix(
        1: n - 1,
        nrow = n,
        ncol = n,
        byrow = TRUE
      ) -
        (1: n - 1))
      rho ^ exponent
    }
    
    R <- ar1_cor(n_pred, alpha)
    return (R)
  }
  
  # # independence
  return (diag(rep(1, times = n_pred)))
}

tmp_R <- build_corr_X(n_pred = 10, corstr = "block diagonal",
                      blocksize = c(2, 4, 2, 1, 2), alpha = 0.25)

tmp_R < -build_corr_X(n_pred = 3, corstr = "block diagonal",
                      blocksize = c(2, 2), alpha = 0.25)

