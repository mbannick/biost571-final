rm(list = ls())

library(nlme)
library(gee)
library(MASS)
library(geepack)
library(extraDistr)
library(expm)
source("fit-gee.R")

set.seed(10)

################################################################
#########  very simple, varing number of observations
################################################################

## three covariates (including intercept)
## beta_0=8,beta_1=-4,beta_2=-2
m<-200
sim_normal_ex <- function(n_sub = m,
                          min_obs = 1,
                          max_obs = 6) {

  # overall setting
  a = 8
  b1 <- (-4)
  b2 <- (-4)
  sigma <- 1
  alpha <- 0.5

  id_vec<-vector()
  a_vec<-vector()
  b1_vec<-vector()
  b2_vec<-vector()
  measure_idx_vec<-vector()
  err_vec<-vector()

  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(
      1:n - 1,
      nrow = n,
      ncol = n,
      byrow = TRUE
    ) -
      (1:n - 1))
    rho ^ exponent
  }

  # individual subject
  for (id in 1:n_sub){

    ni <- rdunif(1, min_obs, max_obs)

    # id index
    id_vec<- c(id_vec,
               rep(id, each = ni))

    # measurement index
    measure_idx_vec <- c(measure_idx_vec,
                         c(1:ni))

    # intercept
    a_vec <- c(a_vec,
               rep(1, times = ni))

    # b1
    b1_vec <- c(b1_vec ,
                sample(c(0,1),ni,replace = TRUE))

    # b2
    b2_vec <- c(b2_vec,
                runif(ni))

    # exchangeable
    #R <- ar1_cor(ni, 0)*sigma
    R <- matrix(alpha,nrow=ni,ncol=ni)
    diag(R) <- 1
    R <- R* sigma

    err <- mvrnorm(n = 1,
                   mu = rep(0, times = ni),
                   Sigma = R)
    err <- as.vector(t(err))
    err_vec<-c(err_vec, err)


  }

  # outcome
  outcome <- a*a_vec + b1 * b1_vec + b2 * b2_vec + err_vec

  long_dat <- data.frame(
    id = id_vec,
    x0 = a_vec,
    x1 = b1_vec,
    x2 = b2_vec,
    outcome = outcome,
    measure_idx=measure_idx_vec,
    error = err_vec
  )


  return(long_dat)
}


data <- sim_normal_ex()

fit <- gee(outcome ~ x1 + x2,
           id = id,
           data = data,
           corstr = "exchangeable")

# FIT GEE WITH OUR NEW FUNCTION

X <- as.matrix(data[,c('x0','x1','x2')])
Y <- data$outcome

source("fit-gee.R")

fit.gee(y=Y, X=X, id=data$id, link='identity', corstr='exchangeable', verbose=T)

fit.gee(y=Y, X=X, id=data$id, link='identity', corstr='exchangeable', verbose=T,
        constrain=TRUE, rho=100)

fit.gee(y=Y, X=X, id=data$id, link='identity', corstr='exchangeable', verbose=T,
        constrain=TRUE, rho=100, gamma=1, lambda=5, tol=5e-4)

fit.gee(y=Y, X=X, id=data$id, link='identity', corstr='exchangeable', verbose=T,
        constrain=TRUE, rho=100, gamma=2, lambda=10, tol=5e-4)
