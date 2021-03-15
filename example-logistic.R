rm(list = ls())

library(nlme)
library(gee)
library(MASS)
library(geepack)
library(extraDistr)
library(expm)
library(data.table)
library(SimCorMultRes)
library(bindata)
source("fit-gee.R")

################################################################
#########  very simple, varing number of observations
################################################################
## three covariates (including intercept)
## beta_0=8,beta_1=-4,beta_2=-2
m<-100
n<-2
alpha <- 0.1

# Get the marginal means within the two groups
x1 <- rnorm(m*n)
x2 <- rnorm(m*n)

# Covariance matrix within subject i
cor.matrix <- toeplitz(c(1, rep(alpha, n-1)))

# Simulate binary outcomes
set.seed(100)
simulated_binary_dataset <- rbin(clsize=n,
                                 intercepts=c(-2), betas=c(0.5, 1.5),
                                 xformula=~x1+x2, cor.matrix=cor.matrix, link="logit")

data <- simulated_binary_dataset$simdata %>% data.table
outcome <- data$y

fit<-gee(outcome ~ x1 + x2,
         id = id,
         data = data,
         family= binomial,
         corstr = "independence")

coef(fit)

# FIT GEE WITH OUR NEW FUNCTION

covars <- as.matrix(data[, c("x1", "x2")])
X <- cbind(rep(1, m*n), covars)
Y <- data$y

# PLAIN
fit.gee(y=Y, X=X, id=data$id, link='logit', corstr='independence', verbose=T,
        constrain=FALSE, rho=100)

# CONSTRAINED
fit.gee(y=Y, X=X, id=data$id, link='logit', corstr='independence', verbose=T,
        constrain=TRUE, rho=100, tol=1e-15)

# LASSO
fit.gee(y=Y, X=X, id=data$id, link='logit', corstr='independence', verbose=T,
        constrain=FALSE, rho=100, tol=1e-15,
        gamma=1, lambda=100)

# RIDGE
fit.gee(y=Y, X=X, id=data$id, link='logit', corstr='independence', verbose=T,
        constrain=FALSE, rho=100, tol=1e-15,
        gamma=2, lambda=10)

# CONSTRAINED LASSO
fit.gee(y=Y, X=X, id=data$id, link='logit', corstr='independence', verbose=T,
        constrain=TRUE, rho=100, tol=1e-15,
        gamma=1, lambda=0.5)

# CONSTRAINED RIDGE
fit.gee(y=Y, X=X, id=data$id, link='logit', corstr='independence', verbose=T,
        constrain=TRUE, rho=100, tol=1e-15,
        gamma=2, lambda=10)
