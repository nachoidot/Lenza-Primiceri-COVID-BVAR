# irf.R
# Impulse response function computation
# Exact port from bvarIrfs.m (25 lines)

#' Compute impulse response functions
#'
#' Computes IRFs using Cholesky ordering to orthogonalize shocks.
#' Based on VAR coefficients and residual covariance matrix.
#'
#' @param beta VAR coefficient matrix (k x n) where k = 1 + n*lags
#' @param sigma Residual covariance matrix (n x n)
#' @param nshock Position of shock variable (1 to n)
#' @param hmax Maximum horizon for IRFs
#'
#' @return Matrix (hmax x n) of impulse responses
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After estimation
#' irfs <- bvar_irfs(result$postmax$betahat,
#'                   result$postmax$sigmahat,
#'                   nshock = 1,
#'                   hmax = 60)
#' }
#'
#' @note Direct port from MATLAB bvarIrfs.m (25 lines)
bvar_irfs <- function(beta, sigma, nshock, hmax) {

  k <- nrow(beta)
  n <- ncol(beta)
  lags <- (k - 1) / n

  ## IRFs computation (lines 14-24 in MATLAB)
  cholVCM <- t(chol(sigma))  # Lower triangular Cholesky factor
  Y <- matrix(0, nrow = lags + hmax, ncol = n)
  in_idx <- lags

  vecshock <- rep(0, n)
  vecshock[nshock] <- 1

  for (tau in 1:hmax) {
    xT <- as.vector(t(Y[(in_idx + tau - 1):(in_idx + tau - lags), , drop = FALSE]))
    xT <- matrix(xT, nrow = 1)

    Y[in_idx + tau, ] <- xT %*% beta[2:k, ] + (tau == 1) * t(cholVCM %*% vecshock)
  }

  irf <- Y[(in_idx + 1):nrow(Y), , drop = FALSE]

  return(irf)
}
