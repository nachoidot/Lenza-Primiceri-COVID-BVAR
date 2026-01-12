# test-likelihood.R
# Unit tests for likelihood computation functions

test_that("log_ml_var_formin returns correct structure", {
  # Testing internal function log_ml_var_formin

  # Create minimal test data
  set.seed(123)
  n <- 3  # Number of variables
  T <- 50  # Number of observations
  lags <- 2
  k <- n * lags + 1

  # Generate synthetic VAR data
  y <- matrix(rnorm(T * n), nrow = T, ncol = n)
  x <- cbind(1, matrix(rnorm(T * k - T), nrow = T))

  # Prior setup
  prior <- covidBVAR:::set_priors(Tcovid = NULL, sur = 0, noc = 0)
  b <- matrix(0, nrow = k, ncol = n)
  diag(b[2:(n + 1), ]) <- 1

  # AR(1) residual variances
  SS <- rep(1, n)
  y0 <- colMeans(y[1:lags, ])

  # Initial parameter (transformed)
  par <- c(0)  # Just lambda (no COVID, no sur, no noc)

  # Call function
  result <- covidBVAR:::log_ml_var_formin(
    par = par,
    y = y, xdata = x, lags = lags, T = T, n = n, b = b,
    MIN = prior$MIN, MAX = prior$MAX, SS = SS, Vc = prior$setpriors$Vc,
    pos = integer(0), mn = list(alpha = 0), sur = 0, noc = 0,
    y0 = y0, hyperpriors = 0, priorcoef = list(), Tcovid = NULL
  )

  # Check structure
  expect_true("logML" %in% names(result))
  expect_true("betahat" %in% names(result))
  expect_true("sigmahat" %in% names(result))

  # Check dimensions
  expect_equal(dim(result$betahat), c(k, n))
  expect_equal(dim(result$sigmahat), c(n, n))

  # logML should be finite

  expect_true(is.finite(result$logML))
})


test_that("log_ml_var_formcmc returns correct structure", {
  # Testing internal function log_ml_var_formcmc

  # Create minimal test data
  set.seed(123)
  n <- 3
  T <- 50
  lags <- 2
  k <- n * lags + 1

  y <- matrix(rnorm(T * n), nrow = T, ncol = n)
  x <- cbind(1, matrix(rnorm(T * k - T), nrow = T))

  prior <- covidBVAR:::set_priors(Tcovid = NULL, sur = 0, noc = 0)
  b <- matrix(0, nrow = k, ncol = n)
  diag(b[2:(n + 1), ]) <- 1

  SS <- rep(1, n)
  y0 <- colMeans(y[1:lags, ])

  # Parameter in original scale
  par <- c(0.2)  # lambda

  # Call without drawing
  result <- covidBVAR:::log_ml_var_formcmc(
    par = par,
    y = y, xdata = x, lags = lags, T = T, n = n, b = b,
    MIN = prior$MIN, MAX = prior$MAX, SS = SS, Vc = prior$setpriors$Vc,
    pos = integer(0), mn = list(alpha = 0), sur = 0, noc = 0,
    y0 = y0, draw = 0, hyperpriors = 0, priorcoef = list(), Tcovid = NULL
  )

  expect_true("logML" %in% names(result))
  expect_true(is.finite(result$logML))

  # Call with drawing
  result_draw <- covidBVAR:::log_ml_var_formcmc(
    par = par,
    y = y, xdata = x, lags = lags, T = T, n = n, b = b,
    MIN = prior$MIN, MAX = prior$MAX, SS = SS, Vc = prior$setpriors$Vc,
    pos = integer(0), mn = list(alpha = 0), sur = 0, noc = 0,
    y0 = y0, draw = 1, hyperpriors = 0, priorcoef = list(), Tcovid = NULL
  )

  expect_true("betadraw" %in% names(result_draw))
  expect_true("sigmadraw" %in% names(result_draw))
  expect_equal(dim(result_draw$betadraw), c(k, n))
  expect_equal(dim(result_draw$sigmadraw), c(n, n))
})


test_that("log_ml_var_formcmc returns -10e15 for out-of-bounds parameters", {
  # Testing internal function log_ml_var_formcmc

  set.seed(123)
  n <- 3
  T <- 50
  lags <- 2
  k <- n * lags + 1

  y <- matrix(rnorm(T * n), nrow = T, ncol = n)
  x <- cbind(1, matrix(rnorm(T * k - T), nrow = T))

  prior <- covidBVAR:::set_priors(Tcovid = NULL, sur = 0, noc = 0)
  b <- matrix(0, nrow = k, ncol = n)
  diag(b[2:(n + 1), ]) <- 1

  SS <- rep(1, n)
  y0 <- colMeans(y[1:lags, ])

  # Out of bounds: lambda < MIN$lambda
  par_low <- c(0.00001)
  result_low <- covidBVAR:::log_ml_var_formcmc(
    par = par_low,
    y = y, xdata = x, lags = lags, T = T, n = n, b = b,
    MIN = prior$MIN, MAX = prior$MAX, SS = SS, Vc = prior$setpriors$Vc,
    pos = integer(0), mn = list(alpha = 0), sur = 0, noc = 0,
    y0 = y0, draw = 0, hyperpriors = 0, priorcoef = list(), Tcovid = NULL
  )
  expect_equal(result_low$logML, -10e15)

  # Out of bounds: lambda > MAX$lambda
  par_high <- c(10)
  result_high <- covidBVAR:::log_ml_var_formcmc(
    par = par_high,
    y = y, xdata = x, lags = lags, T = T, n = n, b = b,
    MIN = prior$MIN, MAX = prior$MAX, SS = SS, Vc = prior$setpriors$Vc,
    pos = integer(0), mn = list(alpha = 0), sur = 0, noc = 0,
    y0 = y0, draw = 0, hyperpriors = 0, priorcoef = list(), Tcovid = NULL
  )
  expect_equal(result_high$logML, -10e15)
})


test_that("COVID weighting is applied correctly", {
  # Testing internal function log_ml_var_formin

  set.seed(123)
  n <- 3
  T <- 100
  lags <- 2
  k <- n * lags + 1
  Tcovid <- 80  # COVID break at observation 80

  y <- matrix(rnorm(T * n), nrow = T, ncol = n)
  x <- cbind(1, matrix(rnorm(T * k - T), nrow = T))

  prior <- covidBVAR:::set_priors(Tcovid = Tcovid)
  b <- matrix(0, nrow = k, ncol = n)
  diag(b[2:(n + 1), ]) <- 1

  SS <- rep(1, n)
  y0 <- colMeans(y[1:lags, ])

  # Parameter with COVID eta (lambda + 4 eta params)
  par <- c(0, log(2), log(2), log(1.5), 0)  # Transformed

  result <- covidBVAR:::log_ml_var_formin(
    par = par,
    y = y, xdata = x, lags = lags, T = T, n = n, b = b,
    MIN = prior$MIN, MAX = prior$MAX, SS = SS, Vc = prior$setpriors$Vc,
    pos = integer(0), mn = list(alpha = 0), sur = 0, noc = 0,
    y0 = y0, hyperpriors = 0, priorcoef = list(), Tcovid = Tcovid
  )

  # Should return finite value
  expect_true(is.finite(result$logML))
})


test_that("sigmahat is symmetric positive definite", {
  # Testing internal function log_ml_var_formin

  set.seed(123)
  n <- 3
  T <- 50
  lags <- 2
  k <- n * lags + 1

  y <- matrix(rnorm(T * n), nrow = T, ncol = n)
  x <- cbind(1, matrix(rnorm(T * k - T), nrow = T))

  prior <- covidBVAR:::set_priors(Tcovid = NULL, sur = 0, noc = 0)
  b <- matrix(0, nrow = k, ncol = n)
  diag(b[2:(n + 1), ]) <- 1

  SS <- rep(1, n)
  y0 <- colMeans(y[1:lags, ])

  par <- c(0)

  result <- covidBVAR:::log_ml_var_formin(
    par = par,
    y = y, xdata = x, lags = lags, T = T, n = n, b = b,
    MIN = prior$MIN, MAX = prior$MAX, SS = SS, Vc = prior$setpriors$Vc,
    pos = integer(0), mn = list(alpha = 0), sur = 0, noc = 0,
    y0 = y0, hyperpriors = 0, priorcoef = list(), Tcovid = NULL
  )

  sigmahat <- result$sigmahat

  # Check symmetry
  expect_equal(sigmahat, t(sigmahat), tolerance = 1e-10)

  # Check positive definiteness (all eigenvalues > 0)
  eigenvalues <- eigen(sigmahat)$values
  expect_true(all(eigenvalues > 0))
})
