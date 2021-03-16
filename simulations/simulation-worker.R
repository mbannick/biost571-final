######################
# Simulation Script
######################

rm(list=ls())

library(data.table)
library(magrittr)
library(R.utils)

setwd("~/repos/biost571-final")
source("fit-gee.R")
source("data-continuous.R")

# GET COMMAND-LINE ARGUMENTS ----------------------

args <- commandArgs(trailingOnly=TRUE, asValues=TRUE,
                    defaults=list(

                      # RANDOM SEED
                      sim=1,

                      # ADD IN ARGUMENTS FOR SIMULATION
                      n_sub=50,
                      min_obs=1,
                      max_obs=6,
                      n_pred=40,
                      sigma_X=1,
                      sigma_Y=2,
                      corr_X='independence',
                      beta_mean=2,
                      alpha_obs=0.0,

                      # THESE ARE ARGUMENTS FOR MODEL FITTING

                      fit_link='identity',
                      corstr='independence',
                      rho=NULL,
                      gamma=NULL,
                      lambda=NULL,
                      tol=1e-12,
                      maxiter=500,
                      verbose=FALSE,
                      directory="."
                    ))

print(args)

# CONVERT COMMAND-LINE ARGUMENTS ------------------

if(is.null(args$rho)){
  args$rho <- NA
  args$constrain <- FALSE
} else {
  args$rho <- as.numeric(args$rho)
  args$constrain <- TRUE
}

if(is.null(args$gamma)){
  args$gamma <- NA
} else {
  args$gamma <- as.numeric(args$gamma)
}

if(is.null(args$lambda)){
  args$lambda <- NA
} else {
  args$lambda <- as.numeric(args$lambda)
}

if(args$corr_X == "block") args$corr_X <- "block diagonal"

# SIMULATE DATA -----------------------------------

# Set a reproducible seed
set.seed(args$sim)

bl.size <- 10
alpha_X <- 0.5

# Simulate data
corr_X <- build_corr_X(n_pred=args$n_pred,
                       alpha=alpha_X,
                       corstr=args$corr_X,
                       blocksize=rep(bl.size, args$n_pred/bl.size))

# Simulate some betas
b.N <- rnorm(n=args$n_pred/4, mean=args$beta_mean)
b.N2 <- -b.N
b.Y <- rep(0, args$n_pred/2)
betas <- c(b.N, b.N2, b.Y)
betas <- sample(betas)

if(args$fit_link == "identity"){

  data <- sim_normal(
    n_sub=args$n_sub,
    min_obs=args$min_obs, max_obs=args$max_obs,
    n_pred=args$n_pred,
    sigma_X=args$sigma_X,
    sigma_Y=args$sigma_Y,
    corr_X=corr_X,
    alpha_obs=args$alpha_obs,
    beta=betas
  )

  y <- data$Y
  X <- cbind(rep(1, length(y)), data$X)
  id <- data$id

}

# FIT MODEL ---------------------------------------

# Fit the GEE model with arguments to simulated data
fit <- fit.gee(y=y, X=X, id=id,
               link=args$fit_link, corstr=args$corstr,
               constrain=args$constrain,
               rho=args$rho, gamma=args$gamma,
               lambda=args$lambda,
               tol=args$tol, maxiter=args$maxiter,
               verbose=args$verbose)

# FIT MODEL W/ CROSS-VALIDATION -------------------

# Sample 70% of the individuals
unique.ids <- unique(id)
iss <- sample(unique.ids, size=0.7*length(unique.ids))
indices.is <- id %in% iss
indices.os <- !id %in% iss

id.is <- id[indices.is]
y.is <- y[indices.is]
X.is <- X[indices.is,]

fit.oos <- fit.gee(y=y.is, X=X.is, id=id.is,
                   link=args$fit_link, corstr=args$corstr,
                   constrain=args$constrain,
                   rho=args$rho, gamma=args$gamma,
                   lambda=args$lambda,
                   tol=args$tol, maxiter=args$maxiter,
                   verbose=args$verbose)

y.os <- y[indices.os]
X.os <- X[indices.os, ]

# COMPUTE EMPIRICAL PERFORMANCE METRICS -----------

true.beta <- c(0, betas)

# Compute either MSE or classification rate
if(args$fit_link == "identity"){

  residual.is <- y - X %*% fit$beta
  mse.is <- mean(residual.is**2)**0.5

  residual.os <- y.os - X.os %*% fit.oos$beta
  mse.os <- mean(residual.os**2)**0.5

} else {

  nu.is <- exp(X %*% fit$beta)/(1 + exp(X %*% fit$beta))
  class.is <- nu.is >= 0.5
  mse.is <- sum(class.is == y) / length(y)

  nu.os <- exp(X.os %*% fit.oos$beta)/(1 + exp(X.os %*% fit.oos$beta))
  class.os <- nu.os >= 0.5
  mse.os <- sum(class.os == y.os) / length(y.os)

}

# Compute the percent of correctly zeroed out predictors

beta0.true <- true.beta <= 1e-15
beta0.fit <- fit$beta <= 1e-15

pct.0 <- mean(beta0.true == beta0.fit)

# SAVE RESULTS ------------------------------------

# Format the results from the fit
df.results <- data.table(

  # REPRODUCIBLE SEED
  sim_id=args$sim,

  # DATA SIMULATION SETTINGS
  # ADD IN ARGUMENTS FOR SIMULATION
  n_sub=args$n_sub,
  min_obs=args$min_obs,
  max_obs=args$max_obs,
  n_pred=args$n_pred,
  sigma_X=args$sigma_X,
  sigma_Y=args$sigma_Y,
  corr_X=args$corr_X,
  beta_mean=args$beta_mean,
  alpha_obs=args$alpha_obs,
  snr=c(data$SNR),

  # MODEL SETTINGS
  fit_link=args$link,
  corstr=args$corstr,
  constrain=args$constrain,
  rho=args$rho,
  gamma=args$gamma,
  lambda=args$lambda,
  tol=args$tol,
  maxiter=args$maxiter,

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
  pct0=pct.0
)

# Create a unique filename for the simulation
directory <- args$directory
args$directory <- NULL
filename <- do.call(paste, args)
filename <- gsub(" ", "_", filename)

# Save data
write.csv(df.results, paste0(directory, filename, ".csv"))
