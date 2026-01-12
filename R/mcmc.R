# mcmc.R
# MCMC sampler using Metropolis-Hastings algorithm
# Exact port from bvarGLP_covid.m lines 194-412

#' Run MCMC for BVAR posterior
#'
#' Implements Metropolis-Hastings algorithm for sampling from posterior
#' distribution of hyperparameters and VAR coefficients.
#'
#' This is extracted from bvarGLP_covid.m (lines 194-412) for modularity.
#' Algorithm and all numerical details are identical to MATLAB.
#'
#' @param y Data matrix (after dummy obs if applicable)
#' @param xdata Regressor matrix (after dummy obs if applicable)
#' @param lags Number of lags
#' @param T Number of observations
#' @param n Number of variables
#' @param b Prior mean on coefficients
#' @param MIN List of minimum bounds
#' @param MAX List of maximum bounds
#' @param SS AR(1) residual variances
#' @param Vc Prior variance on constant
#' @param pos Variables in first differences
#' @param mn List with $alpha
#' @param sur SUR prior flag
#' @param noc NOC prior flag
#' @param y0 Initial condition
#' @param hyperpriors Hyperprior flag
#' @param priorcoef Hyperprior coefficients
#' @param Tcovid COVID break position
#' @param postmode Posterior mode (from optimization)
#' @param xh Optimal parameters (transformed)
#' @param M Total MCMC draws
#' @param N Burn-in draws
#' @param const Proposal scaling constant
#' @param MCMCfcast Generate forecasts flag
#' @param MCMCstorecoeff Store coefficients flag
#' @param hz Forecast horizons
#'
#' @return List with MCMC results
#'
#' @note Direct port from MATLAB bvarGLP_covid.m lines 194-412
#'
#' @noRd
run_mcmc <- function(y, xdata, lags, T, n, b, MIN, MAX, SS, Vc, pos,
                     mn, sur, noc, y0, hyperpriors, priorcoef, Tcovid,
                     postmode, xh, M, N, const, MCMCfcast, MCMCstorecoeff, hz) {

  # Get number of hyperparameters
  k <- 1 + n * lags

  # Recover posterior mode (in original scale)
  if (!is.null(Tcovid)) {
    modeeta <- postmode$eta
  } else {
    modeeta <- numeric(0)
  }

  if (mn$alpha == 1) {
    modealpha <- postmode$alpha
  } else {
    modealpha <- numeric(0)
  }

  if (sur == 1) {
    modetheta <- postmode$theta
  } else {
    modetheta <- numeric(0)
  }

  if (noc == 1) {
    modemiu <- postmode$miu
  } else {
    modemiu <- numeric(0)
  }

  postmode_vec <- c(postmode$lambda, modeeta, modetheta, modemiu, modealpha)

  ## Compute Hessian at posterior mode (using numDeriv package)
  cat('Computing Hessian for MCMC proposal...\n')

  fun <- function(par) {
    result <- log_ml_var_formcmc(par, y, xdata, lags, T, n, b, MIN, MAX, SS, Vc, pos,
                                  mn, sur, noc, y0, 0, hyperpriors, priorcoef, Tcovid)
    return(result$logML)
  }

  Hess <- numDeriv::hessian(fun, postmode_vec)

  # Regularize Hessian (matching MATLAB lines 261-275)
  # If Hessian is not negative definite, use eigenvalue regularization
  eig_decomp <- eigen(Hess)
  V <- Re(eig_decomp$vectors)
  E <- Re(eig_decomp$values)

  # Invert Hessian with robustness guard
  HH <- tryCatch({
    # Attempt direct inversion
    -solve(Hess)
  }, error = function(e) {
    # If direct inversion fails, use eigenvalue regularization
    # Set a floor for eigenvalues (negative eigenvalues indicate positive definite -Hess)
    # For -Hess to be positive definite, Hess eigenvalues should be negative
    cat('Hessian near-singular, applying eigenvalue regularization...\n')

    # Take absolute value of eigenvalues (as MATLAB does for non-PD matrices)
    E_abs <- abs(E)
    # Set minimum eigenvalue floor
    E_floor <- max(abs(E)) * 1e-10
    E_reg <- pmax(E_abs, E_floor)

    # Reconstruct regularized inverse: HH = -Hess^{-1} = V * diag(1/(-E)) * V'
    # But since we want positive definite HH, use: HH = V * diag(1/E_reg) * V'
    V %*% diag(1 / E_reg) %*% t(V)
  })

  # Ensure HH is symmetric (numerical stability)
  HH <- (HH + t(HH)) / 2

  # Ensure HH is positive definite for proposal distribution
  eig_HH <- eigen(HH)
  E_HH <- Re(eig_HH$values)
  if (any(E_HH <= 0)) {
    cat('Proposal covariance not positive definite, fixing eigenvalues...\n')
    V_HH <- Re(eig_HH$vectors)
    E_HH[E_HH <= 0] <- abs(E_HH[E_HH <= 0]) + 1e-10
    HH <- V_HH %*% diag(E_HH) %*% t(V_HH)
    HH <- (HH + t(HH)) / 2
  }

  # Special handling for limited post-COVID data (lines 276-280 in MATLAB)
  if (!is.null(Tcovid) && T <= Tcovid + 1) {
    HessNew <- Hess
    HessNew[4, ] <- 0
    HessNew[, 4] <- 0
    HessNew[4, 4] <- -1

    HH <- tryCatch({
      -solve(HessNew)
    }, error = function(e) {
      cat('Special Hessian inversion failed, using regularized version...\n')
      # Use eigenvalue decomposition for robustness
      eig_new <- eigen(HessNew)
      V_new <- Re(eig_new$vectors)
      E_new <- Re(eig_new$values)
      E_new_abs <- abs(E_new)
      E_new_floor <- max(abs(E_new)) * 1e-10
      E_new_reg <- pmax(E_new_abs, E_new_floor)
      V_new %*% diag(1 / E_new_reg) %*% t(V_new)
    })

    HH[4, 4] <- HH[2, 2]
  }

  ## Storage matrices
  P <- matrix(0, nrow = M, ncol = length(xh))
  LOGML <- numeric(M)

  if (MCMCstorecoeff == 1) {
    beta_store <- array(0, dim = c(k, n, M - N))
    sigma_store <- array(0, dim = c(n, n, M - N))
  } else {
    beta_store <- NULL
    sigma_store <- NULL
  }

  if (MCMCfcast == 1) {
    Dforecast <- array(0, dim = c(length(hz), n, M - N))
  } else {
    Dforecast <- NULL
  }

  ## Starting value for Metropolis algorithm
  logMLold <- -10e15
  attempt <- 0
  max_attempts <- 100

  while (logMLold == -10e15 && attempt < max_attempts) {
    attempt <- attempt + 1
    P[1, ] <- MASS::mvrnorm(1, postmode_vec, HH * const^2)

    result <- log_ml_var_formcmc(P[1, ], y, xdata, lags, T, n, b, MIN, MAX, SS, Vc, pos,
                                  mn, sur, noc, y0,
                                  max(MCMCfcast, MCMCstorecoeff),
                                  hyperpriors, priorcoef, Tcovid)
    logMLold <- result$logML
    betadrawold <- result$betadraw
    sigmadrawold <- result$sigmadraw
  }

  if (logMLold == -10e15) {
    stop("Failed to find valid starting point for MCMC after ", max_attempts, " attempts")
  }

  LOGML[1] <- logMLold

  ## Metropolis-Hastings iterations
  count <- 0

  for (i in 2:M) {
    if (i == 100 * floor(0.01 * i)) {
      cat(sprintf('Now running the %dth mcmc iteration (out of %d)\n', i, M))
    }

    # Propose new value
    P[i, ] <- MASS::mvrnorm(1, P[i - 1, ], HH * const^2)

    # Evaluate proposal
    result_new <- log_ml_var_formcmc(P[i, ], y, xdata, lags, T, n, b, MIN, MAX, SS, Vc, pos,
                                      mn, sur, noc, y0,
                                      max(MCMCfcast, MCMCstorecoeff),
                                      hyperpriors, priorcoef, Tcovid)
    logMLnew <- result_new$logML
    betadrawnew <- result_new$betadraw
    sigmadrawnew <- result_new$sigmadraw

    LOGML[i] <- logMLnew

    ## Metropolis acceptance step
    if (logMLnew > logMLold) {
      # Accept
      logMLold <- logMLnew
      count <- count + 1
    } else {
      if (runif(1) < exp(logMLnew - logMLold)) {
        # Accept with probability
        logMLold <- logMLnew
        count <- count + 1
      } else {
        # Reject
        P[i, ] <- P[i - 1, ]
        LOGML[i] <- logMLold

        # If MCMCfcast or MCMCstorecoeff, redraw beta/sigma with old hyperparams
        if (MCMCfcast == 1 || MCMCstorecoeff == 1) {
          result_old <- log_ml_var_formcmc(P[i, ], y, xdata, lags, T, n, b, MIN, MAX, SS, Vc, pos,
                                            mn, sur, noc, y0,
                                            max(MCMCfcast, MCMCstorecoeff),
                                            hyperpriors, priorcoef, Tcovid)
          betadrawnew <- result_old$betadraw
          sigmadrawnew <- result_old$sigmadraw
        }
      }
    }

    ## Store draws (after burn-in)
    if (i > N) {
      if (MCMCstorecoeff == 1) {
        beta_store[, , i - N] <- betadrawnew
        sigma_store[, , i - N] <- sigmadrawnew
      }

      # Generate forecasts if requested
      if (MCMCfcast == 1) {
        Y <- rbind(y, matrix(0, nrow = max(hz), ncol = n))

        for (tau in 1:max(hz)) {
          xT <- c(1, as.vector(t(Y[(T + tau - 1):(T + tau - lags), , drop = FALSE])))
          xT <- matrix(xT, nrow = 1)

          # COVID volatility scaling (lines 349-357 in MATLAB)
          if (!is.null(Tcovid)) {
            if (T == Tcovid) {
              scaling <- P[i, 3] * (tau == 1) +
                         (1 + (P[i, 4] - 1) * P[i, 5]^(tau - 2)) * (tau >= 2)
            } else if (T > Tcovid) {
              scaling <- 1 + (P[i, 4] - 1) * P[i, 5]^(T - Tcovid + tau - 2)
            }
          } else {
            scaling <- 1
          }

          errors <- MASS::mvrnorm(1, rep(0, n), sigmadrawnew) * scaling
          Y[T + tau, ] <- xT %*% betadrawnew + errors
        }

        Dforecast[, , i - N] <- Y[T + hz, , drop = FALSE]
      }
    }
  }

  ## Extract hyperparameter draws (after burn-in)
  lambda_draws <- P[(N + 1):M, 1]

  if (!is.null(Tcovid)) {
    ncp <- 4
    eta_draws <- P[(N + 1):M, 2:(ncp + 1)]

    if (sur == 1) {
      theta_draws <- P[(N + 1):M, ncp + 2]
      if (noc == 1) {
        miu_draws <- P[(N + 1):M, ncp + 3]
      } else {
        miu_draws <- NULL
      }
    } else {
      theta_draws <- NULL
      if (noc == 1) {
        miu_draws <- P[(N + 1):M, ncp + 2]
      } else {
        miu_draws <- NULL
      }
    }
  } else {
    eta_draws <- NULL
    if (sur == 1) {
      theta_draws <- P[(N + 1):M, 2]
      if (noc == 1) {
        miu_draws <- P[(N + 1):M, 3]
      } else {
        miu_draws <- NULL
      }
    } else {
      theta_draws <- NULL
      if (noc == 1) {
        miu_draws <- P[(N + 1):M, 2]
      } else {
        miu_draws <- NULL
      }
    }
  }

  if (mn$alpha == 1) {
    alpha_draws <- P[(N + 1):M, ncol(P)]
  } else {
    alpha_draws <- NULL
  }

  ## Compute acceptance rate
  ACCrate <- mean(lambda_draws[2:length(lambda_draws)] != lambda_draws[1:(length(lambda_draws) - 1)])

  cat(sprintf('MCMC acceptance rate: %.3f\n', ACCrate))

  ## Return structure
  mcmc_result <- list(
    lambda = lambda_draws,
    eta = eta_draws,
    theta = theta_draws,
    miu = miu_draws,
    alpha = alpha_draws,
    LOGML = LOGML[(N + 1):M],
    beta = beta_store,
    sigma = sigma_store,
    Dforecast = Dforecast,
    ACCrate = ACCrate,
    HH = HH
  )

  return(mcmc_result)
}
