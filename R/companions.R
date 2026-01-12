# companions.R
# Companion matrix formation for VAR state-space representation
# Exact port of FormCompanionMatrices nested function (Baseline_May2021.m lines 122-148)

#' Form VAR Companion Matrices for State-Space Representation
#'
#' Constructs the matrices for the VAR companion form with time-varying
#' COVID volatility. Used for conditional forecasting via Kalman filter.
#'
#' @param betadraw VAR coefficient matrix (k x n) where k = n*lags + 1
#' @param G Cholesky factor of residual covariance (n x n, lower triangular)
#' @param etapar COVID volatility parameters (length 4): [eta1, eta2, eta3, eta4]
#' @param tstar Time period of first COVID shock (March 2020 in forecast period)
#' @param n Number of variables
#' @param lags Number of VAR lags
#' @param TTfcst Total forecast horizon
#'
#' @return List with companion form matrices:
#' \describe{
#'   \item{varc}{Observation constants (n x TTfcst), all zeros}
#'   \item{varZ}{Observation matrices (n x ns x TTfcst) where ns = n*lags}
#'   \item{varG}{Observation noise matrices (n x n x TTfcst), all zeros}
#'   \item{varC}{State constant (ns x 1), includes VAR intercept}
#'   \item{varT}{State transition matrix (ns x ns), companion form}
#'   \item{varH}{State noise matrices (ns x n x TTfcst), time-varying for COVID}
#' }
#'
#' @details
#' State-space representation:
#' \deqn{y_t = c + Z * s_t + G * eps_t}
#' \deqn{s_t = C + T * s_{t-1} + H_t * eta_t}
#'
#' where s_t = [y_t', y_{t-1}', ..., y_{t-lags+1}']' (companion form)
#'
#' COVID volatility scaling in H_t:
#' - t < tstar: H_t = [G', 0, ..., 0]'
#' - t = tstar: H_t = [eta1 * G', 0, ..., 0]'
#' - t = tstar + 1: H_t = [eta2 * G', 0, ..., 0]'
#' - t = tstar + 2: H_t = [eta3 * G', 0, ..., 0]'
#' - t > tstar + 2: H_t = [(1 + (eta3-1)*eta4^(t-tstar-2)) * G', 0, ..., 0]'
#'
#' @note Direct port from MATLAB Baseline_May2021.m lines 122-148
#'
#' @noRd
form_companion_matrices <- function(betadraw, G, etapar, tstar, n, lags, TTfcst) {

  # Observation equation matrices
  # y_t = c + Z * s_t + G * eps_t
  varc <- matrix(0, nrow = n, ncol = TTfcst)  # No observation constant

  # Z matrix: selects first n elements of state (current y)
  # Z = [I_n, 0, ..., 0]  (n x n*lags)
  varZ <- matrix(0, nrow = n, ncol = n * lags)
  varZ[, 1:n] <- diag(n)
  # Replicate for all time periods
  varZ <- array(rep(varZ, TTfcst), dim = c(n, n * lags, TTfcst))

  # No observation noise
  varG <- array(0, dim = c(n, n, TTfcst))

  # State equation matrices
  # s_t = C + T * s_{t-1} + H_t * eta_t
  # where s_t = [y_t', y_{t-1}', ..., y_{t-lags+1}']'

  B <- betadraw

  # State constant: [intercept', 0, ..., 0]'
  varC <- matrix(0, nrow = n * lags, ncol = 1)
  varC[1:n, ] <- B[1, ]  # VAR intercept

  # State transition matrix (companion form)
  # T = [B2:end'          ]  (n x n*lags, VAR coefficients without intercept)
  #     [I_{n*(lags-1)}, 0]  (n*(lags-1) x n*lags, identity for lagged values)
  varT <- matrix(0, nrow = n * lags, ncol = n * lags)
  varT[1:n, ] <- t(B[2:nrow(B), ])  # VAR lag coefficients
  if (lags > 1) {
    varT[(n + 1):(n * lags), 1:(n * (lags - 1))] <- diag(n * (lags - 1))
  }

  # State noise matrices (time-varying for COVID)
  # H_t = [scaling * G']
  #       [0           ]
  #       [...         ]
  #       [0           ]  (n*lags x n)
  varH <- array(0, dim = c(n * lags, n, TTfcst))

  for (t in 1:TTfcst) {
    if (t < tstar) {
      # Pre-COVID period: no scaling
      varH[1:n, , t] <- G
    } else if (t == tstar) {
      # March 2020: eta1 scaling
      varH[1:n, , t] <- G * etapar[1]
    } else if (t == tstar + 1) {
      # April 2020: eta2 scaling
      varH[1:n, , t] <- G * etapar[2]
    } else if (t == tstar + 2) {
      # May 2020: eta3 scaling
      varH[1:n, , t] <- G * etapar[3]
    } else {
      # Post-May 2020: geometric decay to 1
      # scaling = 1 + (eta3 - 1) * eta4^(t - tstar - 2)
      scaling <- 1 + (etapar[3] - 1) * etapar[4]^(t - tstar - 2)
      varH[1:n, , t] <- G * scaling
    }
  }

  return(list(
    varc = varc,
    varZ = varZ,
    varG = varG,
    varC = varC,
    varT = varT,
    varH = varH
  ))
}
