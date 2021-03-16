######################
# Simulation Script
######################

rm(list=ls())

library(data.table)
library(magrittr)
library(R.utils)

setwd("~/repos/biost571-final")
source("fit-gee.R")
source("sim_dat_with_logistic.R")

# GET COMMAND-LINE ARGUMENTS ----------------------

args <- commandArgs(trailingOnly=TRUE)
sim <- args[1]
# sim <- 1

set.seed(sim)

# SIMULATE DATA -----------------------------------

# FIXED SIMULATION PARAMETERS
sigma_X <- 1
sigma_Y <- 2
beta_mean <- 2
rho <- 100
tol <- 1e-2
maxiter <- 250

# SIMULATION PARAMETERS
n_sub.list <- c(100)
min_obs.list <- c(2)
str_X.list <- c("independence", "block")
fit_link.list <- c("identity", "logit")
alpha_obs.list <- c(0.4)
n_pred.list <- c(80)
gamma.list <- c(NA, 1, 2)
lambda.list <- c(0, 0.5, 1, 1.5)
constrain.list <- c(F, T)

run.simulation <- function(n_sub, min_obs, str_X,
                           fit_link, alpha_obs,
                           n_pred, gamma, lambda, constrain){

  bl.size <- 10
  alpha_X <- 0.5

  # Simulate data
  corr_X <- build_corr_X(n_pred=n_pred,
                         alpha=alpha_X,
                         corstr=str_X,
                         blocksize=rep(bl.size, n_pred/bl.size))

  # Simulate some betas
  b.N <- rnorm(n=n_pred/4, mean=beta_mean)
  b.N2 <- -b.N
  b.Y <- rep(0, n_pred/2)
  betas <- c(b.N, b.N2, b.Y)
  betas <- sample(betas)

  data <- sim_norm_logit(
    n_sub=n_sub,
    min_obs=min_obs, max_obs=max_obs,
    n_pred=n_pred,
    sigma_X=sigma_X,
    sigma_Y=sigma_Y,
    corr_X=corr_X,
    alpha_obs=alpha_obs,
    beta=betas,
    link=fit_link
  )

  y <- data$Y
  X <- cbind(rep(1, length(y)), data$X)
  id <- data$id

  # FIT MODEL W/ CROSS-VALIDATION -------------------

  # Sample 70% of the individuals
  unique.ids <- unique(id)
  iss <- sample(unique.ids, size=0.9*length(unique.ids))
  indices.is <- id %in% iss
  indices.os <- !id %in% iss

  id.is <- id[indices.is]
  y.is <- y[indices.is]
  X.is <- X[indices.is,]

  fit <- fit.gee(y=y.is, X=X.is, id=id.is,
                     link=fit_link, corstr=corstr,
                     constrain=constrain,
                     rho=rho, gamma=gamma,
                     lambda=lambda,
                     tol=tol, maxiter=maxiter,
                     verbose=FALSE)

  y.os <- y[indices.os]
  X.os <- X[indices.os, ]

  # COMPUTE EMPIRICAL PERFORMANCE METRICS -----------

  true.beta <- c(0, betas)

  # Compute either MSE or classification rate
  if(fit_link == "identity"){

    residual.is <- y.is - X.is %*% fit$beta
    mse.is <- mean(residual.is**2)**0.5

    residual.os <- y.os - X.os %*% fit$beta
    mse.os <- mean(residual.os**2)**0.5

  } else {

    nu.is <- exp(X.is %*% fit$beta)/(1 + exp(X.is %*% fit$beta))
    class.is <- nu.is >= 0.5
    mse.is <- sum(class.is == y.is) / length(y.is)

    nu.os <- exp(X.os %*% fit$beta)/(1 + exp(X.os %*% fit$beta))
    class.os <- nu.os >= 0.5
    mse.os <- sum(class.os == y.os) / length(y.os)

  }

  # Compute the percent of correctly zeroed out predictors

  beta0.true <- true.beta <= 1e-15
  beta0.fit <- fit$beta <= 1e-15

  pct.0 <- mean(beta0.true == beta0.fit)

  const.sum <- sum(fit$beta)

  # SAVE RESULTS ------------------------------------

  # Format the results from the fit
  results <- data.table(

    # REPRODUCIBLE SEED
    sim_id=sim,

    # DATA SIMULATION SETTINGS
    # ADD IN ARGUMENTS FOR SIMULATION
    n_sub=n_sub,
    min_obs=min_obs,
    max_obs=max_obs,
    n_pred=n_pred,
    sigma_X=sigma_X,
    sigma_Y=sigma_Y,
    str_X=str_X,
    beta_mean=beta_mean,
    alpha_obs=alpha_obs,
    snr=c(data$SNR),

    # MODEL SETTINGS
    fit_link=fit_link,
    corstr=corstr,
    constrain=constrain,
    rho=rho,
    gamma=gamma,
    lambda=lambda,
    tol=tol,
    maxiter=maxiter,

    # MODEL RESULTS

    # Whether the model converged
    success=fit$success,
    # What the observed tolerance was
    error=fit$error,
    # What was the estimate of phi
    phi=fit$phi,
    # What was the estimated alpha
    alpha=fit$alpha,

    # EMPIRICAL PERFORMANCE

    mse_is=mse.is,
    mse_os=mse.os,
    pct0=pct.0,
    sum=const.sum
  )

  return(results)
}

df <- data.table()

# LOOP THROUGH SIMULATION SETTINGS
for(n_sub in n_sub.list){
  for(min_obs in min_obs.list){
    max_obs <- min_obs + 5
    for(str_X in str_X.list){
      for(fit_link in fit_link.list){
        if(fit_link == "identity") corstr <- "exchangeable"
        if(fit_link == "logit") corstr <- "independence"
        for(alpha_obs in alpha_obs.list){
          for(n_pred in n_pred.list){
            for(gamma in gamma.list){
              for(lambda in lambda.list){
                for(con in constrain.list){
                  print(n_sub)
                  print(min_obs)
                  print(str_X)
                  print(fit_link)
                  print(alpha_obs)
                  print(n_pred)
                  print(gamma)
                  print(lambda)
                  print(con)
                  if(is.na(gamma) & lambda != 0){
                    next
                  }
                  if(!is.na(gamma) & lambda == 0){
                    next
                  }
                  df.result <- run.simulation(n_sub=n_sub, min_obs=min_obs,
                                              str_X=str_X,
                                              fit_link=fit_link,
                                              alpha_obs=alpha_obs,
                                              n_pred=n_pred,
                                              gamma=gamma, lambda=lambda,
                                              constrain=con)
                  df <- rbind(df, df.result)
                }
              }
            }
          }
        }
      }
    }
  }
}

# Create a unique filename for the simulation
directory <- args$directory
args$directory <- NULL
filename <- paste0("simulation_", sim, ".csv")

# Save data
write.csv(df, paste0(directory, filename))
