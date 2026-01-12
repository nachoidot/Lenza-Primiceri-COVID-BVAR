# utils.R
# Utility functions ported from MATLAB code
# All functions are exact translations preserving MATLAB behavior

#' Create lagged matrix
#'
#' Creates a matrix or vector of lagged values, exactly replicating MATLAB's lag.m function
#'
#' @param x Input matrix or vector (nobs x k)
#' @param n Order of lag (default = 1)
#' @param v Initial values for padding (default = 0)
#'
#' @return Matrix (or vector) of lags (nobs x k)
#'
#' @details
#' This function exactly replicates the MATLAB lag.m function from the original codebase.
#' If n <= 0, returns NULL (equivalent to MATLAB's []).
#' The first n rows are filled with value v, and the remaining rows contain
#' the original data shifted down by n rows (last n rows are removed).
#'
#' @note Direct port from MATLAB lag.m by James P. LeSage
#'
#' @noRd
lag <- function(x, n = 1, v = 0) {
  # Handle matrix or vector input
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }

  # Get number of columns
  nc <- ncol(x)

  # Handle edge cases
  if (n < 1) {
    return(NULL)  # Equivalent to MATLAB's []
  }

  # Create initial padding with value v
  zt <- matrix(v, nrow = n, ncol = nc)

  # Combine padding with trimmed x
  # trimr(x, 0, n) removes last n rows
  z <- rbind(zt, trimr(x, 0, n))

  return(z)
}


#' Trim rows from matrix
#'
#' Returns a matrix (or vector) x stripped of the specified rows.
#' Exact port of MATLAB trimr.m function.
#'
#' @param x Input matrix (or vector) (n x k)
#' @param n1 First n1 rows to strip
#' @param n2 Last n2 rows to strip
#'
#' @return z = x[(n1+1):(n-n2), ]
#'
#' @details
#' Modeled after Gauss trimr function. Removes n1 rows from the beginning
#' and n2 rows from the end of matrix x.
#'
#' @note Direct port from MATLAB trimr.m by James P. LeSage
#'
#' @noRd
trimr <- function(x, n1, n2) {
  # Handle vector input
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }

  n <- nrow(x)

  # Check if trying to trim too much
  if ((n1 + n2) >= n) {
    stop("Attempting to trim too much in trimr")
  }

  h1 <- n1 + 1
  h2 <- n - n2

  z <- x[h1:h2, , drop = FALSE]

  return(z)
}


#' Return number of columns
#'
#' Returns the number of columns in a matrix x.
#' Exact port of MATLAB cols.m function.
#'
#' @param x Input matrix
#'
#' @return Number of columns in x
#'
#' @note Direct port from MATLAB cols.m by James P. LeSage
#'
#' @noRd
cols <- function(x) {
  # Handle vector input
  if (is.vector(x)) {
    return(1)
  }
  return(ncol(x))
}


#' OLS estimation
#'
#' Performs OLS regression exactly as in MATLAB ols1.m
#'
#' @param y Dependent variable vector (nobs x 1)
#' @param x Independent variables matrix (nobs x nvar)
#'
#' @return List with OLS results:
#' \describe{
#'   \item{nobs}{Number of observations}
#'   \item{nvar}{Number of variables}
#'   \item{bhatols}{OLS coefficient estimates}
#'   \item{yhatols}{Fitted values}
#'   \item{resols}{Residuals}
#'   \item{sig2hatols}{Residual variance (unbiased estimate)}
#'   \item{sigbhatols}{Variance-covariance matrix of coefficients}
#'   \item{XX}{X'X matrix}
#'   \item{R2}{R-squared}
#' }
#'
#' @details
#' This function exactly replicates MATLAB's ols1.m function.
#' Note: sig2hatols uses (nobs - nvar) as denominator for unbiased estimate.
#'
#' @note Direct port from MATLAB ols1.m
#'
#' @noRd
ols1 <- function(y, x) {
  # Convert to matrix if needed
  if (is.vector(y)) {
    y <- matrix(y, ncol = 1)
  }
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }

  nobs <- nrow(x)
  nvar <- ncol(x)
  nobsy <- nrow(y)

  # Check dimensions
  if (nobs != nobsy) {
    stop("x and y must have same # obs in ols")
  }

  # OLS estimation - exactly as MATLAB
  XX <- t(x) %*% x
  bhatols <- solve(XX) %*% (t(x) %*% y)
  yhatols <- x %*% bhatols
  resols <- y - yhatols
  sig2hatols <- as.numeric((t(resols) %*% resols) / (nobs - nvar))
  sigbhatols <- sig2hatols * solve(XX)
  R2 <- var(as.vector(yhatols)) / var(as.vector(y))

  # Return structure matching MATLAB
  r <- list(
    nobs = nobs,
    nvar = nvar,
    bhatols = bhatols,
    yhatols = yhatols,
    resols = resols,
    sig2hatols = sig2hatols,
    sigbhatols = sigbhatols,
    XX = XX,
    R2 = R2
  )

  return(r)
}


#' Gamma distribution coefficients
#'
#' Computes Gamma(k, theta) shape and scale parameters given mode and standard deviation.
#' Exact port of MATLAB GammaCoef.m function.
#'
#' @param mode Mode of the Gamma distribution
#' @param sd Standard deviation of the Gamma distribution
#' @param plotit Flag for plotting (default = 0, not implemented in R version)
#'
#' @return List with:
#' \describe{
#'   \item{k}{Shape parameter}
#'   \item{theta}{Scale parameter}
#' }
#'
#' @details
#' Uses the formulas:
#' k = (2 + mode^2/sd^2 + sqrt((4 + mode^2/sd^2) * mode^2/sd^2)) / 2
#' theta = sqrt(sd^2 / k)
#'
#' @note Direct port from MATLAB GammaCoef.m
#'
#' @noRd
gamma_coef <- function(mode, sd, plotit = 0) {
  # Exact MATLAB formula
  k <- (2 + mode^2/sd^2 + sqrt((4 + mode^2/sd^2) * mode^2/sd^2)) / 2
  theta <- sqrt(sd^2 / k)

  r <- list(
    k = k,
    theta = theta
  )

  # Note: plotit functionality not implemented (not needed for estimation)
  if (plotit == 1) {
    warning("Plotting functionality not implemented in R version")
  }

  return(r)
}
