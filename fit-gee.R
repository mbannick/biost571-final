
#' Fit a GEE model with sum 0 constraint and possibly a LASSO
#' or ridge penalty.
#'
#' @param y Outcome vector
#' @param X Design matrix
#' @param id Vector of ids, should align with design matrix
#' and outcome vector
#' @param link Link function, one of 'identity' or 'logit'
#' @param corstr Correlation structure, one of 'independence' or
#' 'exchangeable'
#' @param constrain Whether or not to impose a sum to 0 constraint
#' on the regression coefficients
#' @param rho The step size for the lagrange multiplier update
#' if using the sum to 0 constraint. Ignored if constrain = FALSE.
#' @param gamma An integer specifying the order of the penalty,
#' 1 corresponds to LASSO and >=2 corresponds to bridge penalty (2 = Ridge)
#' @param lambda The hyper-parameter for penalization. Ignored if is.na(gamma)
#' @param tol Tolerance for the regression coefficients
#' @param maxiter Maximum number of iterations
#' @param verbose Whether or not want information printed at each iteration
fit.gee <- function(y, X, id, link='identity', corstr='independence',
                    constrain=FALSE, rho=NA, gamma=NA, lambda=NA,
                    tol=1e-5, maxiter=500, verbose=FALSE, variance=FALSE){

  # DATA ATTRIBUTES -------------------------------------

  # Number of observations
  n_obs <- length(y)

  # Number of parameters
  p_df <- dim(X)[2]

  # Number of individuals
  m <- length(unique(id))

  ids <- unique(id)

  # CORRELATION STRUCTURE -------------------------------

  indep <- corstr == "independence"
  if(indep){
    print("Using an independence working correlation.")
  } else {
    print("Using an exchangeable working correlation.")
  }

  # LINK FUNCTION + FAMILY ------------------------------

  normal <- link == "identity"
  if(normal){
    print("Using an identity link, Gaussian outcomes")
    fam <- gaussian(link="identity")
  } else {
    print("Using a logit link, binary outcomes")
    fam <- binomial(link="logit")
  }

  if(!normal & !indep) stop("Exchangeable logistic not implemented.")

  # PENALTIES -------------------------------------------

  if(!is.na(gamma) & is.na(lambda)) stop("Need to provide a lambda parameter.")

  lasso <- FALSE
  bridge <- FALSE
  if(!is.na(gamma)){
    lasso <- gamma == 1
    bridge <- gamma >= 2
  }

  if(lasso) print(paste0("Implementing LASSO penalty with lambda = ", lambda))
  if(bridge) print(paste0("Implementing bridge penalty with lambda = ", lambda))

  # CONSTRAINTS -----------------------------------------
  if(constrain){
    print("Implementing sum to 0 constraint.")
    A <- matrix(rep(1, p_df), nrow = 1, ncol = p_df)
    if(is.na(rho)) stop("Need to provide a rho step size
                        for Lagrange multiplier in order to use constraints.")
  }

  # INITIALIZATION --------------------------------------
  # Initialize the parameter vector with a glm
  mod.init <- glm(y ~ 0 + X, family=fam)
  beta.init <- coef(mod.init)
  cat("Beta: ", beta.init, "\n")

  # Initialization for phi
  # (This is not used for binary outcomes, for binary phi = 1)
  res.init <- y - X %*% beta.init
  if(normal){
    phi.init <- sum(res.init**2)/(n_obs - p_df)
  } else {
    phi.init <- 1
  }

  # Initialization for alpha
  # (This is not used for independence)
  if(!indep){
    alpha <- 0
    resi_sum <- 0
    n_alpha <- 0
  } else {
    alpha <- NA
  }

  # Initialization for lagrange multiplier
  if(constrain){
    eta <- 1
  } else {
    eta <- NA
  }


  # ALGORITHM -------------------------------------------

  error <- 1
  beta <- beta.init
  phi <- phi.init

  iter <- 1
  success <- FALSE

  while(error > tol){
    if(iter == maxiter){
      break
    }
    if(!verbose) cat(".")

    # Set the score and hessian matrices
    score <- matrix(0, nrow=p_df)
    hess <- matrix(0, nrow=p_df, ncol=p_df)

    for(i in ids){
      # SUBSET THE DATA ---------------------------------
      indices <- id == i
      X_i <- X[indices, ]
      if(sum(indices) == 1) X_i <- t(as.matrix(X_i))
      y_i <- y[indices]
      ni <- length(y_i)

      # Get the linear predictor
      nu_i <- X_i %*% beta

      # BASIC SCORE AND HESSIAN ------------------------
      if(normal){

        # RESIDUAL
        res_i <- y_i - nu_i

        if(indep){

          # ** INDEPENDENT NORMAL **
          # ------------------------

          score_i <- (-t(X_i)) %*% res_i
          hess_i <- t(X_i) %*% X_i

        } else {

          # ** EXCHANGEABLE NORMAL **
          # -------------------------

          # Calculate the variance
          # matrix with correlations
          R_i <- matrix(alpha, nrow=ni, ncol=ni)
          diag(R_i) <- 1
          Vm_i <- diag(ni) * phi
          V_i <- sqrtm(Vm_i) %*% R_i %*% sqrtm(Vm_i)

          score_i <- (-t(X_i)) %*% solve(V_i) %*% res_i
          hess_i <- t(X_i) %*% solve(V_i) %*% X_i

        }
      } else {

        # ** INDEPENDENT LOGISTIC **
        # --------------------------

        # MEAN FUNCTION
        mu_i <- exp(nu_i)/(1 + exp(nu_i))
        mu_i2 <- mu_i * (1 - mu_i)

        # RESIDUAL
        res_i <- y_i - mu_i

        # SCORE
        score_i <- (-t(X_i)) %*% res_i

        # HESSIAN
        if(ni==1){
          hess_i <- t(X_i) %*% X_i * as.numeric(mu_i2)
        } else {
          hess_i <- t(X_i) %*% diag(as.vector(mu_i2)) %*% X_i
        }

      }

      score <- score + score_i
      hess <- hess + hess_i

      # ALPHA: PAIRWISE RESIDUALS ----------------------
      p_resi <- (y_i - X_i %*% beta)/sqrt(phi)

      if(!indep){
        if(ni > 1){
          pairs <- combn(1:ni,2)
          resi_sum <- resi_sum + sum(p_resi[pairs[1,]] * p_resi[pairs[2,]])
          n_alpha <- n_alpha + dim(pairs)[2]
        }
      }

    }

    # LAGRANGE MULTIPLIER: UPDATE SCORE + HESS ---------
    if(constrain){
      score <- score + t(A) %*% eta + rho * t(A) %*% A %*% beta
      hess <- hess + rho * t(A) %*% A
    }

    # RIDGE PENALTY: UPDATE SCORE + HESS ---------------
    if(bridge){
      score <- score + lambda * gamma * sign(beta) * abs(beta)**(gamma - 1)
      hess <- hess + diag(c(lambda * gamma * (gamma - 1) * abs(beta)**(gamma - 2)))
    }

    # NEWTON RAPHSON: UPDATE BETA ----------------------
    beta.new <- beta - solve(hess) %*% score
    if(any(is.nan(beta.new))) break("NaNs in betas. Diverged.")

    # LASSO: SOFT THRESHOLDING ON THE BETAS ------------
    if(lasso){
      beta.new <- ifelse(abs(beta.new) < lambda, 0, beta.new - sign(beta.new) * lambda)
    }

    # METHOD OF MOMENTS: UPDATE ALPHA ------------------
    if(!indep){
      alpha <- resi_sum/(n_alpha - p_df)
    }

    # METHOD OF MOMENTS: UPDATE PHI --------------------
    if(normal){
      resi <- y - X %*% beta.new
      phi <- sum(resi**2)/(n_obs - p_df)
    }

    # LAGRANGE MULTIPLIER: UPDATE ETA ------------------
    if(constrain){
      eta <- eta + rho * A %*% beta.new
    }

    if(verbose){
      cat("Beta: ", beta.new, "\n")
      if(normal) cat("Phi: ", phi, "\n")
      if(!indep) cat("Alpha: ", alpha, "\n")
      if(constrain) cat("Sum to 0: ", sum(beta.new), "\n")
      if(constrain) cat("Eta: ", eta, "\n")
      if(lasso) cat("L1 Norm: ", sum(abs(beta.new)), "\n")
      if(bridge) cat("Bridge: ", sum(abs(beta.new)**gamma), "\n")
    }

    error <- sqrt(sum((beta.new - beta)**2))
    if(verbose) cat("Error: ", error, "\n")
    beta <- beta.new
    iter <- iter + 1

    if(error <= tol) success <- TRUE

  }

  # SANDWICH VARIANCE ESTIMATION ----------------------


  if(variance){
    cat("Computing sandwich variance", "\n")
    model_based <- matrix(data=0, nrow=p_df, ncol=p_df)
    empirical <- matrix(data=0, nrow=p_df, ncol=p_df)

    for(i in ids){

      # SUBSET THE DATA ---------------------------------
      indices <- id == i
      X_i <- X[indices, ]
      if(sum(indices) == 1) X_i <- t(as.matrix(X_i))
      y_i <- y[indices]
      ni <- length(y_i)

      # Get the linear predictor
      nu_i <- X_i %*% beta

      if(normal){

        # RESIDUAL
        res_i <- y_i - nu_i

        if(indep){

          # ** INDEPENDENT NORMAL **
          # ------------------------

          # Calculate model-based and empirical variance
          model_based <- model_based + t(X_i) %*% X_i
          empirical <- empirical + t(X_i) %*% res_i %*% t(res_i) %*% X_i

        } else {

          # ** EXCHANGEABLE NORMAL **
          # -------------------------

          # Calculate the variance
          # matrix with correlations
          R_i <- matrix(alpha, nrow=ni, ncol=ni)
          diag(R_i) <- 1
          Vm_i <- diag(ni) * phi
          V_i <- sqrtm(Vm_i) %*% R_i %*% sqrtm(Vm_i)

          # inverse to save computational time
          V_ii <- solve(V_i)

          # Calculate model-based and empirical variance
          model_based <- model_based + t(X_i) %*% V_ii %*% X_i
          empirical <- empirical + t(X_i) %*% V_ii %*% res_i %*% t(res_i) %*% V_ii %*% X_i

        }
      } else {

        # ** INDEPENDENT LOGISTIC **
        # --------------------------

        # MEAN FUNCTION
        mu_i <- exp(nu_i)/(1 + exp(nu_i))
        mu_i2 <- mu_i * (1 - mu_i)

        # RESIDUAL
        res_i <- y_i - mu_i

        if(ni==1){
          Dt <- as.numeric(mu_i2) * X_i
        } else {
          Dt <- diag(as.vector(mu_i2)) %*% X_i
        }

        # Calculate model-based and empirical variance
        model_based <- model_based + t(X_i) %*% Dt
        empirical <- empirical + t(X_i) %*% res_i %*% t(res_i) %*% X_i

      }
    }

    covar <- solve(model_based) %*% empirical %*% solve(model_based)
  } else {
    covar <- NA
  }

  print("Exiting.")
  if(success){
    cat("Converged with parameters", beta, "\n")
  } else {
    cat("Failed to converge after", maxiter, "iterations!", "\n")
    cat("Achieved error: ", error, "\n")
    cat("Last iteration:", beta, "\n")
  }

  return(list(
    beta=beta,
    alpha=alpha,
    eta=eta,
    phi=phi,
    success=success,
    error=error,
    vcov=covar
  ))

}
