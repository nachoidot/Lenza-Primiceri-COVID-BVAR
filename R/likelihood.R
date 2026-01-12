# likelihood.R
# Log-posterior computation functions
# Exact ports of logMLVAR_formin_covid.m and logMLVAR_formcmc_covid.m

#' Compute log marginal likelihood for optimization
#'
#' Computes the log-posterior (or logML if hyperpriors=0), the posterior mode
#' of the coefficients and the covariance matrix of the residuals of the BVAR
#' of Giannone, Lenza and Primiceri (2015), augmented with a change in
#' volatility at the time of COVID (March 2020).
#'
#' This version is for use during optimization (hyperparameters are logit-transformed).
#'
#' @param par Transformed hyperparameters (logit scale)
#' @param y Data matrix (T x n)
#' @param xdata Regressor matrix (T x k)
#' @param lags Number of lags
#' @param T Number of observations
#' @param n Number of variables
#' @param b Prior mean on VAR coefficients (k x n)
#' @param MIN List of minimum bounds for hyperparameters
#' @param MAX List of maximum bounds for hyperparameters
#' @param SS AR(1) residual variances (n x 1)
#' @param Vc Prior variance on constant
#' @param pos Variables entering in first differences
#' @param mn List with $alpha element (0 or 1)
#' @param sur Sum-of-coefficients prior flag (0 or 1)
#' @param noc No-cointegration prior flag (0 or 1)
#' @param y0 Initial condition (mean of first lags observations)
#' @param hyperpriors Flag for hyperpriors (0 or 1)
#' @param priorcoef Hyperprior coefficients
#' @param Tcovid Position of March 2020 (or NULL)
#'
#' @return List with:
#' \describe{
#'   \item{logML}{Negative log marginal likelihood (for minimization)}
#'   \item{betahat}{Posterior mode of VAR coefficients}
#'   \item{sigmahat}{Posterior mode of covariance matrix}
#' }
#'
#' @note Direct port from MATLAB logMLVAR_formin_covid.m (177 lines)
#'
#' @noRd
log_ml_var_formin <- function(par, y, xdata, lags, T, n, b, MIN, MAX, SS, Vc, pos,
                               mn, sur, noc, y0, hyperpriors, priorcoef, Tcovid) {

  ## Hyperparameters - inverse logit transformation
  lambda <- MIN$lambda + (MAX$lambda - MIN$lambda) / (1 + exp(-par[1]))

  d <- n + 2
  psi <- SS * (d - n - 1)

  # Initialize theta, miu, eta
  theta <- MIN$theta
  miu <- MIN$miu
  eta <- MIN$eta

  if (is.null(Tcovid)) {
    # No COVID period
    if (sur == 1) {
      theta <- MIN$theta + (MAX$theta - MIN$theta) / (1 + exp(-par[2]))
      if (noc == 1) {
        miu <- MIN$miu + (MAX$miu - MIN$miu) / (1 + exp(-par[3]))
      }
    } else if (sur == 0) {
      if (noc == 1) {
        miu <- MIN$miu + (MAX$miu - MIN$miu) / (1 + exp(-par[2]))
      }
    }
  } else {
    # COVID period included
    ncp <- 4  # number of COVID parameters

    eta <- MIN$eta + (MAX$eta - MIN$eta) / (1 + exp(-par[2:(ncp+1)]))

    # Construct inverse weights (volatility scaling)
    invweights <- rep(1, T)
    invweights[Tcovid] <- eta[1]
    invweights[Tcovid + 1] <- eta[2]

    if (T > Tcovid + 1) {
      idx <- (Tcovid + 2):T
      invweights[idx] <- 1 + (eta[3] - 1) * eta[4]^(0:(T - Tcovid - 2))
    }

    # Rescale data by inverse weights
    y <- diag(1 / invweights) %*% y
    xdata <- diag(1 / invweights) %*% xdata

    # Extract remaining hyperparameters
    if (sur == 1) {
      theta <- MIN$theta + (MAX$theta - MIN$theta) / (1 + exp(-par[ncp + 2]))
      if (noc == 1) {
        miu <- MIN$miu + (MAX$miu - MIN$miu) / (1 + exp(-par[ncp + 3]))
      }
    } else if (sur == 0) {
      if (noc == 1) {
        miu <- MIN$miu + (MAX$miu - MIN$miu) / (1 + exp(-par[ncp + 2]))
      }
    }
  }

  # Alpha parameter (lag decay)
  if (mn$alpha == 0) {
    alpha <- 2
  } else if (mn$alpha == 1) {
    alpha <- MIN$alpha + (MAX$alpha - MIN$alpha) / (1 + exp(-par[length(par)]))
  }

  ## Setting up the priors
  k <- 1 + n * lags
  omega <- rep(0, k)
  omega[1] <- Vc

  for (i in 1:lags) {
    idx <- (1 + (i - 1) * n + 1):(1 + i * n)
    omega[idx] <- (d - n - 1) * (lambda^2) * (1 / (i^alpha)) / psi
  }

  # Prior scale matrix for the covariance of the shocks
  PSI <- diag(psi)

  # Dummy observations
  Td <- 0
  xdsur <- NULL
  ydsur <- NULL
  xdnoc <- NULL
  ydnoc <- NULL

  if (sur == 1) {
    xdsur <- matrix(c(1/theta, (1/theta) * rep(y0, lags)), nrow = 1)
    ydsur <- matrix((1/theta) * y0, nrow = 1)

    y <- rbind(y, ydsur)
    xdata <- rbind(xdata, xdsur)
    Td <- 1
  }

  if (noc == 1) {
    ydnoc <- (1/miu) * diag(as.vector(y0))
    if (length(pos) > 0) {
      ydnoc[pos, pos] <- 0
    }

    xdnoc <- cbind(matrix(0, nrow = n, ncol = 1),
                   (1/miu) * matrix(rep(diag(as.vector(y0)), lags), nrow = n))

    y <- rbind(y, ydnoc)
    xdata <- rbind(xdata, xdnoc)
    Td <- Td + n
  }

  T <- T + Td

  ## Output: posterior mode of VAR coefficients
  XtX <- t(xdata) %*% xdata
  betahat <- solve(XtX + diag(1/omega)) %*% (t(xdata) %*% y + diag(1/omega) %*% b)

  # VAR residuals
  epshat <- y - xdata %*% betahat

  # Posterior mode of the covariance matrix
  sigmahat <- (t(epshat) %*% epshat + PSI +
               t(betahat - b) %*% diag(1/omega) %*% (betahat - b)) / (T + d + n + 1)

  ## Log marginal likelihood
  aaa <- diag(sqrt(omega)) %*% XtX %*% diag(sqrt(omega))
  bbb <- diag(1/sqrt(psi)) %*% (t(epshat) %*% epshat +
               t(betahat - b) %*% diag(1/omega) %*% (betahat - b)) %*% diag(1/sqrt(psi))

  # Eigenvalue regularization (exact MATLAB procedure)
  eigaaa <- Re(eigen(aaa, only.values = TRUE)$values)
  eigaaa[eigaaa < 1e-12] <- 0
  eigaaa <- eigaaa + 1

  eigbbb <- Re(eigen(bbb, only.values = TRUE)$values)
  eigbbb[eigbbb < 1e-12] <- 0
  eigbbb <- eigbbb + 1

  # Log-ML formula (lines 121-122 in MATLAB)
  logML <- -n * T * log(pi) / 2 +
           sum(lgamma((T + d - (0:(n-1))) / 2) - lgamma((d - (0:(n-1))) / 2)) +
           -T * sum(log(psi)) / 2 +
           -n * sum(log(eigaaa)) / 2 +
           -(T + d) * sum(log(eigbbb)) / 2

  ## Normalizing constant (if sur or noc)
  if (sur == 1 || noc == 1) {
    yd <- rbind(ydsur, ydnoc)
    xd <- rbind(xdsur, xdnoc)

    # Prior mode of VAR coefficients
    betahatd <- b  # This is the case for our priors

    # VAR residuals at prior mode
    epshatd <- yd - xd %*% betahatd

    XdtXd <- t(xd) %*% xd
    aaa <- diag(sqrt(omega)) %*% XdtXd %*% diag(sqrt(omega))
    bbb <- diag(1/sqrt(psi)) %*% (t(epshatd) %*% epshatd +
                 t(betahatd - b) %*% diag(1/omega) %*% (betahatd - b)) %*% diag(1/sqrt(psi))

    eigaaa <- Re(eigen(aaa, only.values = TRUE)$values)
    eigaaa[eigaaa < 1e-12] <- 0
    eigaaa <- eigaaa + 1

    eigbbb <- Re(eigen(bbb, only.values = TRUE)$values)
    eigbbb[eigbbb < 1e-12] <- 0
    eigbbb <- eigbbb + 1

    # Normalizing constant
    norm <- -n * Td * log(pi) / 2 +
            sum(lgamma((Td + d - (0:(n-1))) / 2) - lgamma((d - (0:(n-1))) / 2)) +
            -Td * sum(log(psi)) / 2 +
            -n * sum(log(eigaaa)) / 2 +
            -(T + d) * sum(log(eigbbb)) / 2

    logML <- logML - norm
  }

  ## Log determinant to account for re-weighting (COVID adjustment)
  if (!is.null(Tcovid)) {
    logML <- logML - n * sum(log(invweights))
  }

  ## Add hyperprior contribution
  if (hyperpriors == 1) {
    logML <- logML + log_gamma_pdf(lambda, priorcoef$lambda$k, priorcoef$lambda$theta)

    if (sur == 1) {
      logML <- logML + log_gamma_pdf(theta, priorcoef$theta$k, priorcoef$theta$theta)
    }
    if (noc == 1) {
      logML <- logML + log_gamma_pdf(miu, priorcoef$miu$k, priorcoef$miu$theta)
    }
    if (!is.null(Tcovid)) {
      logML <- logML - 2 * log(eta[1]) - 2 * log(eta[2]) - 2 * log(eta[3]) +
               log_beta_pdf(eta[4], priorcoef$eta4$alpha, priorcoef$eta4$beta)
    }
  }

  # Return negative logML for minimization
  logML <- -logML

  return(list(logML = logML, betahat = betahat, sigmahat = sigmahat))
}


#' Compute log marginal likelihood for MCMC
#'
#' Computes the log-posterior and draws from the posterior distribution.
#' This version is for MCMC (hyperparameters are in original scale).
#'
#' @param par Hyperparameters in original scale (not transformed)
#' @param ... Same parameters as log_ml_var_formin
#' @param draw Flag: 0 = no draws, 1 = draw from posterior of beta and sigma
#'
#' @return List with:
#' \describe{
#'   \item{logML}{Log marginal likelihood (NOT negated)}
#'   \item{betadraw}{Draw from posterior of beta (if draw=1)}
#'   \item{sigmadraw}{Draw from posterior of sigma (if draw=1)}
#' }
#'
#' @note Direct port from MATLAB logMLVAR_formcmc_covid.m (215 lines)
#'
#' @noRd
log_ml_var_formcmc <- function(par, y, xdata, lags, T, n, b, MIN, MAX, SS, Vc, pos,
                                mn, sur, noc, y0, draw, hyperpriors, priorcoef, Tcovid) {

  ## Hyperparameters (original scale, not transformed)
  lambda <- par[1]

  d <- n + 2
  theta <- MIN$theta
  miu <- MIN$miu
  eta <- MIN$eta
  psi <- SS * (d - n - 1)

  if (is.null(Tcovid)) {
    if (sur == 1) {
      theta <- par[2]
      if (noc == 1) {
        miu <- par[3]
      }
    } else if (sur == 0) {
      if (noc == 1) {
        miu <- par[2]
      }
    }
  } else {
    ncp <- 4
    eta <- par[2:(ncp + 1)]

    # Construct inverse weights
    invweights <- rep(1, T)
    invweights[Tcovid] <- eta[1]
    invweights[Tcovid + 1] <- eta[2]

    if (T > Tcovid + 1) {
      idx <- (Tcovid + 2):T
      invweights[idx] <- 1 + (eta[3] - 1) * eta[4]^(0:(T - Tcovid - 2))
    }

    # Rescale data
    y <- diag(1 / invweights) %*% y
    xdata <- diag(1 / invweights) %*% xdata

    if (sur == 1) {
      theta <- par[ncp + 2]
      if (noc == 1) {
        miu <- par[ncp + 3]
      }
    } else if (sur == 0) {
      if (noc == 1) {
        miu <- par[ncp + 2]
      }
    }
  }

  if (mn$alpha == 0) {
    alpha <- 2
  } else if (mn$alpha == 1) {
    alpha <- par[length(par)]
  }

  ## Return very low posterior if parameters out of bounds
  check_vec <- c(lambda, eta, theta, miu, alpha)
  lower_vec <- c(MIN$lambda, MIN$eta, MIN$theta, MIN$miu, MIN$alpha)
  upper_vec <- c(MAX$lambda, MAX$eta, MAX$theta, MAX$miu, MAX$alpha)

  if (sum(check_vec < lower_vec) > 0 || sum(check_vec > upper_vec) > 0) {
    return(list(logML = -10e15, betadraw = NULL, sigmadraw = NULL))
  }

  ## Priors (same as formin version)
  k <- 1 + n * lags
  omega <- rep(0, k)
  omega[1] <- Vc

  for (i in 1:lags) {
    idx <- (1 + (i - 1) * n + 1):(1 + i * n)
    omega[idx] <- (d - n - 1) * (lambda^2) * (1 / (i^alpha)) / psi
  }

  PSI <- diag(psi)

  Td <- 0
  xdsur <- NULL
  ydsur <- NULL
  xdnoc <- NULL
  ydnoc <- NULL

  if (sur == 1) {
    xdsur <- matrix(c(1/theta, (1/theta) * rep(y0, lags)), nrow = 1)
    ydsur <- matrix((1/theta) * y0, nrow = 1)
    y <- rbind(y, ydsur)
    xdata <- rbind(xdata, xdsur)
    Td <- 1
  }

  if (noc == 1) {
    ydnoc <- (1/miu) * diag(as.vector(y0))
    if (length(pos) > 0) {
      ydnoc[pos, pos] <- 0
    }
    xdnoc <- cbind(matrix(0, nrow = n, ncol = 1),
                   (1/miu) * matrix(rep(diag(as.vector(y0)), lags), nrow = n))
    y <- rbind(y, ydnoc)
    xdata <- rbind(xdata, xdnoc)
    Td <- Td + n
  }

  ## Output: posterior mode
  XtX <- t(xdata) %*% xdata
  betahat <- solve(XtX + diag(1/omega)) %*% (t(xdata) %*% y + diag(1/omega) %*% b)
  epshat <- y - xdata %*% betahat

  T <- T + Td

  ## Log-ML computation
  aaa <- diag(sqrt(omega)) %*% XtX %*% diag(sqrt(omega))
  bbb <- diag(1/sqrt(psi)) %*% (t(epshat) %*% epshat +
               t(betahat - b) %*% diag(1/omega) %*% (betahat - b)) %*% diag(1/sqrt(psi))

  eigaaa <- Re(eigen(aaa, only.values = TRUE)$values)
  eigaaa[eigaaa < 1e-12] <- 0
  eigaaa <- eigaaa + 1

  eigbbb <- Re(eigen(bbb, only.values = TRUE)$values)
  eigbbb[eigbbb < 1e-12] <- 0
  eigbbb <- eigbbb + 1

  logML <- -n * T * log(pi) / 2 +
           sum(lgamma((T + d - (0:(n-1))) / 2) - lgamma((d - (0:(n-1))) / 2)) +
           -T * sum(log(psi)) / 2 +
           -n * sum(log(eigaaa)) / 2 +
           -(T + d) * sum(log(eigbbb)) / 2

  if (sur == 1 || noc == 1) {
    yd <- rbind(ydsur, ydnoc)
    xd <- rbind(xdsur, xdnoc)
    betahatd <- b
    epshatd <- yd - xd %*% betahatd

    XdtXd <- t(xd) %*% xd
    aaa <- diag(sqrt(omega)) %*% XdtXd %*% diag(sqrt(omega))
    bbb <- diag(1/sqrt(psi)) %*% (t(epshatd) %*% epshatd +
                 t(betahatd - b) %*% diag(1/omega) %*% (betahatd - b)) %*% diag(1/sqrt(psi))

    eigaaa <- Re(eigen(aaa, only.values = TRUE)$values)
    eigaaa[eigaaa < 1e-12] <- 0
    eigaaa <- eigaaa + 1

    eigbbb <- Re(eigen(bbb, only.values = TRUE)$values)
    eigbbb[eigbbb < 1e-12] <- 0
    eigbbb <- eigbbb + 1

    norm <- -n * Td * log(pi) / 2 +
            sum(lgamma((Td + d - (0:(n-1))) / 2) - lgamma((d - (0:(n-1))) / 2)) +
            -Td * sum(log(psi)) / 2 +
            -n * sum(log(eigaaa)) / 2 +
            -(T + d) * sum(log(eigbbb)) / 2

    logML <- logML - norm
  }

  if (!is.null(Tcovid)) {
    logML <- logML - n * sum(log(invweights))
  }

  if (hyperpriors == 1) {
    logML <- logML + log_gamma_pdf(lambda, priorcoef$lambda$k, priorcoef$lambda$theta)
    if (sur == 1) {
      logML <- logML + log_gamma_pdf(theta, priorcoef$theta$k, priorcoef$theta$theta)
    }
    if (noc == 1) {
      logML <- logML + log_gamma_pdf(miu, priorcoef$miu$k, priorcoef$miu$theta)
    }
    if (!is.null(Tcovid)) {
      logML <- logML - 2 * log(eta[1]) - 2 * log(eta[2]) - 2 * log(eta[3]) +
               log_beta_pdf(eta[4], priorcoef$eta4$alpha, priorcoef$eta4$beta)
    }
  }

  ## Draw from posterior (if requested)
  betadraw <- NULL
  sigmadraw <- NULL

  if (draw == 1) {
    S <- PSI + t(epshat) %*% epshat +
         t(betahat - b) %*% diag(1/omega) %*% (betahat - b)

    # Draw sigma from Inverse Wishart
    eig_decomp <- eigen(S)
    V <- eig_decomp$vectors
    E <- eig_decomp$values
    Sinv <- V %*% diag(1/abs(E)) %*% t(V)

    # Draw T+d vectors from N(0, Sinv), then invert
    eta_draws <- matrix(rnorm(n * (T + d)), nrow = T + d, ncol = n) %*% chol(Sinv)
    sigmadraw <- solve(t(eta_draws) %*% eta_draws)

    # Draw beta from multivariate normal
    cholSIGMA <- cholred(0.5 * (sigmadraw + t(sigmadraw)))
    cholZZinv <- cholred(solve(XtX + diag(1/omega)))
    betadraw <- betahat + t(cholZZinv) %*% matrix(rnorm(nrow(betahat) * ncol(betahat)),
                                                   nrow = nrow(betahat)) %*% cholSIGMA
  }

  return(list(logML = logML, betadraw = betadraw, sigmadraw = sigmadraw))
}


## Helper functions for log densities

#' Log Gamma PDF
#' @noRd
log_gamma_pdf <- function(x, k, theta) {
  (k - 1) * log(x) - x/theta - k * log(theta) - lgamma(k)
}

#' Log Beta PDF
#' @noRd
log_beta_pdf <- function(x, al, bet) {
  (al - 1) * log(x) + (bet - 1) * log(1 - x) - lbeta(al, bet)
}


#' Regularized Cholesky decomposition
#'
#' Performs Cholesky decomposition with regularization for near-singular matrices.
#' Exact port of MATLAB cholred nested function.
#'
#' @param S Symmetric matrix
#'
#' @return Upper triangular Cholesky factor (transposed to match MATLAB)
#'
#' @note Direct port from MATLAB cholred function in logMLVAR_formcmc_covid.m
#'
#' @noRd
cholred <- function(S) {
  # Symmetrize
  S <- 0.5 * (S + t(S))

  # Eigenvalue decomposition
  eig_decomp <- eigen(S)
  v <- eig_decomp$vectors
  d <- Re(eig_decomp$values)

  # Regularization threshold
  scale <- mean(diag(S)) * 1e-12
  J <- (d > scale)

  # Construct Cholesky factor
  C <- matrix(0, nrow = ncol(S), ncol = nrow(S))
  if (sum(J) > 0) {
    C[J, ] <- t(v[, J] %*% diag(sqrt(d[J]), nrow = sum(J)))
  }

  return(t(C))  # Return transposed to match MATLAB convention
}
