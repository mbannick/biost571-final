######################
# Simulation Script
######################

rm(list=ls())

library(data.table)
library(magrittr)
library(R.utils)

setwd("~/repos/biost571-final")
source("fit-gee.R")

# GET COMMAND-LINE ARGUMENTS ----------------------

args <- commandArgs(trailingOnly=TRUE, asValues=TRUE,
                    defaults=list(

                      # RANDOM SEED
                      sim=1,

                      # ADD IN ARGUMENTS FOR SIMULATION

                      # THESE ARE ARGUMENTS FOR MODEL FITTING

                      fit_link='identity',
                      corstr='independence',
                      constrain=FALSE,
                      rho=NULL,
                      gamma=NULL,
                      lambda=NULL,
                      tol=1e-15,
                      maxiter=500,
                      verbose=FALSE,
                      directory="."
                    ))

# CONVERT COMMAND-LINE ARGUMENTS ------------------

if(is.null(args$rho)){
  args$rho <- NA
} else {
  args$rho <- as.numeric(args$rho)
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

# SIMULATE DATA -----------------------------------

# Set a reproducible seed
set.seed(args$sim)

# Simulate data
# ...

# FIT MODEL ---------------------------------------

# Fit the GEE model with arguments to simulated data
fit <- fit.gee(y=y, X=X, id=id,
               link=args$link, corstr=args$corstr,
               constrain=args$constrain,
               rho=args$rho, gamma=args$gamma,
               lambda=args$lambda,
               tol=args$tol, maxiter=args$maxiter,
               verbose=args$verbose)

# COMPUTE EMPIRICAL PERFORMANCE METRICS -----------

# Decide what we want to compute, i.e. bias, coverage, MSE, etc.
# ...

bias <- 0
coverage <- FALSE
mse <- 0

# SAVE RESULTS ------------------------------------

# Format the results from the fit
df.results <- data.table(

  # REPRODUCIBLE SEED
  sim_id=sim,

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

  bias=bias,
  coverage=coverage,
  mse=mse
)

# Create a unique filename for the simulation
filename <- do.call(paste, args)
filename <- gsub(" ", "_", filename)

# Save data
write.csv(df.results, paste0(args$directory, filename, ".csv"))
