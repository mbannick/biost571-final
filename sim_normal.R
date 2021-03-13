rm(list = ls())

library(nlme)
library(gee)
library(MASS)
library(geepack)
library(extraDistr)
library(expm)

#set.seed(10)


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
b2 <- (-2)
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


simp_dat <- sim_normal_ex()


################################################################
#########  very simple, varing number of observations
################################################################

## geefit
fit <-
  geese(
    outcome ~ x1 + x2,
    id = id,
    data = simp_dat,
    sformula = ~ x1,
    corstr = "ar1",
    jack = TRUE,
    j1s = TRUE,
    fij = TRUE
  )
summary(fit)

if(FALSE){
fit<-gee(outcome ~ x1 + x2,
         id = id,
         data = simp_dat,
         corstr = "exchangeable")
}

################################################################
######### solve by hand (using GEE, no constriant)
################################################################

# initialize
beta <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
X <- as.matrix(simp_dat[,c('x0','x1','x2')])
Y <- simp_dat$outcome
phi <- sd(Y)
n_obs <- length(Y)
p_df <- length(beta)
alpha <-0

# iterate
for (q in 1:10){
  cat(".")
  # update phi
  resi<-Y-X %*% beta
  p_resi<-(Y-X %*% beta)/sqrt(phi)
  phi<-sum(resi^2)/(n_obs-p_df)

  # for updating alpha
  resi_sum <- 0
  n_alpha <- 0

  # update GEE
  score <- matrix(0,nrow=length(beta),ncol=1)
  hes_gee <- matrix(0,nrow=length(beta),ncol=length(beta))

  for ( id in 1: m){
    # individual contribution to the score and hessian

    # subject specific data
    one_sub<-simp_dat[simp_dat$id==id,]
    X_i <- as.matrix(one_sub[, c('x0','x1','x2')])
    Y_i <- one_sub$outcome
    ni <- length(Y_i)

    # working correlation
    R_i <- matrix(alpha,nrow=ni,ncol=ni)
    diag(R_i) <- 1
    Vm_i <- diag(ni)*phi
    V_i <- sqrtm(Vm_i) %*% R_i %*% sqrtm(Vm_i)

    # score and hessian
    resi_i <- Y_i - X_i %*% beta
    score_i <- (-t(X_i)) %*% solve(V_i) %*% resi_i
    hes_gee_i <- t(X_i) %*% solve(V_i) %*% X_i

    score <- score+score_i
    hes_gee <- hes_gee+hes_gee_i

    # for updating alpha
    p_resi<-(Y_i-X_i %*% beta)/sqrt(phi)

    if(ni>1){
    pairs <- combn(1:ni,2)
    resi_sum  <- resi_sum+
                sum(p_resi[pairs[1,]]*p_resi[pairs[2,]])
    n_alpha <- n_alpha+ dim(pairs)[2]
    }

  }

  # update beta
  beta <- beta - solve(hes_gee) %*% score

  # update alpha
  alpha <- resi_sum/(n_alpha-p_df)

}


beta
sum(beta)



################################################################
######### solve by hand (using GEE, with sum to zero constriant)
################################################################


# initialize
beta <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
A <- matrix(rep(1, 3), nrow = 1, ncol = 3)
eta <- 1
rho <- 100 ## too small can't converge?

# ADDED BRIDGE PENALTY
# If you want to run non-penalized regression, set lambda = 0
lambda <- 5
gamma <- 2 # corresponding to ridge regression

X <- as.matrix(simp_dat[,c('x0','x1','x2')])
Y <- simp_dat$outcome
n_obs <- length(Y)
p_df <- length(beta)

phi <- sd(Y)
alpha <- 0

# iterate
for (q in 1:30){

  # update phi
  resi<-Y-X %*% beta
  p_resi<-(Y-X %*% beta)/sqrt(phi)
  phi<-sum(resi^2)/(n_obs-p_df)

  # for updating alpha
  resi_sum <- 0
  n_alpha <- 0

  # update gee
  score_gee <- matrix(0,nrow=length(beta),ncol=1)
  hes_gee <- matrix(0,nrow=length(beta),ncol=length(beta))

  for ( id in 1:m){

    # individual contribution to the score and hessian

    # subject specific data
    one_sub<-simp_dat[simp_dat$id==id,]
    X_i <- as.matrix(one_sub[, c('x0','x1','x2')])
    Y_i <- one_sub$outcome
    ni <- length(Y_i)

    # working correlation
    R_i <- matrix(alpha,nrow=ni,ncol=ni)
    diag(R_i) <- 1
    Vm_i <- diag(ni)*phi
    V_i <- sqrtm(Vm_i) %*% R_i %*% sqrtm(Vm_i)

    # score and hessian
    resi_i <- Y_i - X_i %*% beta
    score_i <- (-t(X_i)) %*% solve(V_i) %*% resi_i
    hes_gee_i <- t(X_i) %*% solve(V_i) %*% X_i

    score_gee <- score_gee+score_i
    hes_gee <- hes_gee+hes_gee_i

    # for updating alpha
    p_resi<-(Y_i-X_i %*% beta)/sqrt(phi)

    if(ni>1){
      pairs <- combn(1:ni,2)
      resi_sum  <- resi_sum+
        sum(p_resi[pairs[1,]]*p_resi[pairs[2,]])
      n_alpha <- n_alpha+ dim(pairs)[2]
    }

  }

  score <- score_gee + t(A) %*% eta + rho * t(A) %*% A %*% beta
  # Add in the bridge penalty
  bridge.score <- lambda * gamma * sign(beta) * abs(beta)**(gamma - 1)
  score <- score + bridge.score
  hes <- hes_gee + rho * t(A) %*% A
  bridge.hes <- diag(c(lambda * gamma * (gamma - 1) * abs(beta)**(gamma - 2)))
  # note that we technically have sign(beta) * sign(beta), but this is always 1
  hes <- hes + bridge.hes

  # update beta and eta
  beta <- beta - solve(hes) %*% score
  eta <- eta + rho * A %*% beta

  # update alpha
  alpha <- resi_sum/(n_alpha-p_df)

  print(beta)
  print(sum(beta))
  print(sum(abs(beta)**gamma))
}

beta
sum(beta)
sum(abs(beta))
sum(beta**2)
