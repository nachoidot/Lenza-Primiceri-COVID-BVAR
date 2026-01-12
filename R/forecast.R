# forecast.R
# VAR forecasting functions
# Exact port from bvarFcst.m (22 lines)

#' Compute VAR forecasts
#'
#' Computes forecasts for VAR using specified coefficients.
#' This is unconditional forecasting (no shocks).
#'
#' @param y Historical data matrix (T x n)
#' @param beta VAR coefficient matrix (k x n)
#' @param hz Forecast horizons (vector)
#'
#' @return Matrix (length(hz) x n) of forecasts
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After estimation
#' forecasts <- bvar_forecast(y, result$postmax$betahat, hz = 1:12)
#' }
#'
#' @note Direct port from MATLAB bvarFcst.m (22 lines)
bvar_forecast <- function(y, beta, hz) {

  k <- nrow(beta)
  n <- ncol(beta)
  lags <- (k - 1) / n
  T <- nrow(y)

  ## Compute forecasts (lines 13-20 in MATLAB)
  Y <- rbind(y, matrix(0, nrow = max(hz), ncol = n))

  for (tau in 1:max(hz)) {
    xT <- c(1, as.vector(t(Y[(T + tau - 1):(T + tau - lags), , drop = FALSE])))
    xT <- matrix(xT, nrow = 1)
    Y[T + tau, ] <- xT %*% beta
  }

  forecast <- Y[T + hz, , drop = FALSE]

  return(forecast)
}
