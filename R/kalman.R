# kalman.R
# Kalman filter and disturbance smoother for state-space models
# Exact ports from DisturbanceSmootherVAR.m and kfilter_forDS_VAR.m

#' Kalman filter step
#'
#' Single step of Kalman filter for state-space model.
#' Exact port of kfilter_forDS_VAR.m (17 lines)
#'
#' @param y Observation vector (n x 1)
#' @param c Observation constant (n x 1)
#' @param Z Observation matrix (n x ns)
#' @param G Observation noise matrix (n x n)
#' @param C State constant (ns x 1)
#' @param TT State transition matrix (ns x ns)
#' @param H State noise matrix (ns x ne)
#' @param shat Current state estimate (ns x 1)
#' @param sig Current state covariance (ns x ns)
#'
#' @return List with:
#' \describe{
#'   \item{shatnew}{Updated state estimate}
#'   \item{signew}{Updated state covariance}
#'   \item{v}{Innovation}
#'   \item{k}{Kalman gain}
#'   \item{sigmainv}{Innovation precision}
#' }
#'
#' @note Direct port from MATLAB kfilter_forDS_VAR.m
#'
#' @noRd
kfilter_step <- function(y, c, Z, G, C, TT, H, shat, sig) {
  # Note: TT is called T in MATLAB, renamed to avoid conflict with TRUE

  n <- length(y)

  # Prediction step
  omega <- TT %*% sig %*% t(TT) + H %*% t(H)
  sigmainv <- solve(Z %*% omega %*% t(Z) + G %*% t(G))

  # Kalman gain
  k <- omega %*% t(Z) %*% sigmainv

  # Innovation
  v <- y - c - Z %*% (C + TT %*% shat)

  # Update step
  shatnew <- C + TT %*% shat + k %*% v
  signew <- omega - k %*% Z %*% omega

  return(list(
    shatnew = shatnew,
    signew = signew,
    v = v,
    k = k,
    sigmainv = sigmainv
  ))
}


#' Disturbance smoother for VAR
#'
#' Performs Kalman filtering and disturbance smoothing for state-space model.
#' Can operate in "smoother" or "simulation" mode for conditional forecasting.
#'
#' Exact port of DisturbanceSmootherVAR.m (112 lines)
#'
#' @param y Observation matrix (T x n)
#' @param c Observation constants (n x T)
#' @param Z Observation matrices (n x ns x T)
#' @param G Observation noise matrices (n x n x T)
#' @param C State constant (ns x 1)
#' @param B State transition matrix (ns x ns)
#' @param H State noise matrices (ns x ne x T)
#' @param s00 Initial state (ns x 1)
#' @param P00 Initial state covariance (ns x ns)
#' @param TT Number of time periods
#' @param n Number of observables
#' @param ns Number of states
#' @param ne Number of state noise dimensions
#' @param SS Mode: "smoother" or "simulation"
#'
#' @return List with:
#' \describe{
#'   \item{sdraw}{Smoothed/simulated states (ns x T)}
#'   \item{epsdraw}{Smoothed/simulated shocks (ne x T) or NULL}
#' }
#'
#' @note Direct port from MATLAB DisturbanceSmootherVAR.m (112 lines)
#'
#' @noRd
disturbance_smoother_var <- function(y, c, Z, G, C, B, H, s00, P00, TT, n, ns, ne, SS) {

  # Ensure matrices are 3D arrays
  if (length(dim(Z)) == 2) {
    Z <- array(Z, dim = c(nrow(Z), ncol(Z), TT))
  }
  if (is.matrix(c)) {
    if (ncol(c) == 1) {
      c <- matrix(rep(c, TT), nrow = n, ncol = TT)
    }
  }
  if (length(dim(G)) == 2) {
    G <- array(G, dim = c(nrow(G), ncol(G), TT))
  }
  if (length(dim(H)) == 2) {
    H <- array(H, dim = c(nrow(H), ncol(H), TT))
  }

  # Handle missing values
  ind <- !is.na(y)

  ## Kalman filter (lines 16-36 in MATLAB)
  V <- matrix(0, nrow = n, ncol = TT)
  K <- array(0, dim = c(ns, n, TT))
  HINV <- array(0, dim = c(n, n, TT))
  SHAT <- matrix(0, nrow = ns, ncol = TT)
  SIG <- array(0, dim = c(ns, ns, TT))

  shat <- s00
  sig <- P00

  for (t in 1:TT) {
    y_t <- y[t, ind[t, ]]
    c_t <- c[ind[t, ], t]
    Z_t <- Z[ind[t, ], , t, drop = FALSE]
    G_t <- G[ind[t, ], ind[t, ], t, drop = FALSE]

    kf_result <- kfilter_step(y_t, c_t, Z_t, G_t, C, B, H[, , t], shat, sig)

    shat <- kf_result$shatnew
    sig <- kf_result$signew

    SHAT[, t] <- shat
    SIG[, , t] <- sig
    V[ind[t, ], t] <- kf_result$v
    K[, ind[t, ], t] <- kf_result$k
    HINV[ind[t, ], ind[t, ], t] <- kf_result$sigmainv
  }

  ## Disturbance smoother (lines 38-47 in MATLAB)
  epshat <- matrix(0, nrow = ne, ncol = TT)
  r <- matrix(0, nrow = ns, ncol = 1)

  for (t in TT:1) {
    Z_t <- Z[ind[t, ], , t, drop = FALSE]
    H_t <- H[, , t]
    HINV_t <- HINV[ind[t, ], ind[t, ], t, drop = FALSE]
    V_t <- V[ind[t, ], t]
    K_t <- K[, ind[t, ], t, drop = FALSE]

    epshat[, t] <- t(H_t) %*% t(Z_t) %*% HINV_t %*% V_t +
                   t(H_t) %*% t(diag(ns) - K_t %*% Z_t) %*% r

    r <- t(B) %*% t(Z_t) %*% HINV_t %*% V_t +
         t(B) %*% t(diag(ns) - K_t %*% Z_t) %*% r
  }

  if (SS == "smoother") {
    ## Smoothed states (lines 49-57 in MATLAB)
    sdraw <- matrix(0, nrow = ns, ncol = TT)
    sdraw[, 1] <- C + B %*% s00 + H[, , 1] %*% epshat[, 1]

    for (t in 2:TT) {
      sdraw[, t] <- C + B %*% sdraw[, t - 1] + H[, , t] %*% epshat[, t]
    }

    epsdraw <- NULL

  } else if (SS == "simulation") {
    ## Simulation mode (lines 59-110 in MATLAB)

    # Simulate new shocks
    epsplus <- matrix(rnorm(ne * TT), nrow = ne, ncol = TT)
    splus <- matrix(0, nrow = ns, ncol = TT)
    yplus <- matrix(0, nrow = n, ncol = TT)

    splus[, 1] <- C + B %*% s00 + H[, , 1] %*% epsplus[, 1]
    yplus[ind[1, ], 1] <- Z[ind[1, ], , 1, drop = FALSE] %*% splus[, 1]

    for (t in 2:TT) {
      splus[, t] <- C + B %*% splus[, t - 1] + H[, , t] %*% epsplus[, t]
      yplus[ind[t, ], t] <- Z[ind[t, ], , t, drop = FALSE] %*% splus[, t]
    }

    # Kalman filter on simulated data
    Vplus <- matrix(0, nrow = n, ncol = TT)
    Kplus <- array(0, dim = c(ns, n, TT))
    HINVplus <- array(0, dim = c(n, n, TT))

    shat <- s00
    sig <- P00

    for (t in 1:TT) {
      yplus_t <- yplus[ind[t, ], t]
      c_t <- c[ind[t, ], t]
      Z_t <- Z[ind[t, ], , t, drop = FALSE]
      G_t <- G[ind[t, ], ind[t, ], t, drop = FALSE]

      kf_result <- kfilter_step(yplus_t, c_t, Z_t, G_t, C, B, H[, , t], shat, sig)

      shat <- kf_result$shatnew
      sig <- kf_result$signew
      Vplus[ind[t, ], t] <- kf_result$v
      Kplus[, ind[t, ], t] <- kf_result$k
      HINVplus[ind[t, ], ind[t, ], t] <- kf_result$sigmainv
    }

    # Disturbance smoother on simulated data
    epshatplus <- matrix(0, nrow = ne, ncol = TT)
    r <- matrix(0, nrow = ns, ncol = 1)

    for (t in TT:1) {
      Z_t <- Z[ind[t, ], , t, drop = FALSE]
      H_t <- H[, , t]
      HINVplus_t <- HINVplus[ind[t, ], ind[t, ], t, drop = FALSE]
      Vplus_t <- Vplus[ind[t, ], t]
      Kplus_t <- Kplus[, ind[t, ], t, drop = FALSE]

      epshatplus[, t] <- t(H_t) %*% t(Z_t) %*% HINVplus_t %*% Vplus_t +
                         t(H_t) %*% t(diag(ns) - Kplus_t %*% Z_t) %*% r

      r <- t(B) %*% t(Z_t) %*% HINVplus_t %*% Vplus_t +
           t(B) %*% t(diag(ns) - Kplus_t %*% Z_t) %*% r
    }

    # Conditional simulation (line 104)
    epsdraw <- epshat + epsplus - epshatplus

    # Generate states from conditional shocks
    sdraw <- matrix(0, nrow = ns, ncol = TT)
    sdraw[, 1] <- C + B %*% s00 + H[, , 1] %*% epsdraw[, 1]

    for (t in 2:TT) {
      sdraw[, t] <- C + B %*% sdraw[, t - 1] + H[, , t] %*% epsdraw[, t]
    }
  }

  return(list(sdraw = sdraw, epsdraw = epsdraw))
}
