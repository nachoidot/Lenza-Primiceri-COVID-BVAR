# test-utils.R
# Unit tests for utility functions

test_that("trimr works correctly", {
  # Test case 1: Simple matrix trimming
  x <- matrix(1:10, ncol = 1)
  z <- covidBVAR:::trimr(x, 2, 3)
  expect_equal(z, matrix(3:7, ncol = 1))

  # Test case 2: Matrix with multiple columns
  x <- matrix(1:20, ncol = 2)
  z <- covidBVAR:::trimr(x, 1, 1)
  expect_equal(nrow(z), 8)
  expect_equal(ncol(z), 2)

  # Test case 3: Error when trimming too much
  expect_error(covidBVAR:::trimr(matrix(1:5, ncol = 1), 3, 3))
})


test_that("cols works correctly", {
  # Test case 1: Matrix with multiple columns
  x <- matrix(1:20, ncol = 4)
  expect_equal(covidBVAR:::cols(x), 4)

  # Test case 2: Vector (should return 1)
  x <- c(1, 2, 3, 4, 5)
  expect_equal(covidBVAR:::cols(x), 1)

  # Test case 3: Single column matrix
  x <- matrix(1:10, ncol = 1)
  expect_equal(covidBVAR:::cols(x), 1)
})


test_that("lag works correctly", {
  # Test case 1: Lag 1 with default initial value (0)
  x <- matrix(1:5, ncol = 1)
  z <- covidBVAR:::lag(x, 1)
  expect_equal(z, matrix(c(0, 1, 2, 3, 4), ncol = 1))

  # Test case 2: Lag 2
  x <- matrix(1:5, ncol = 1)
  z <- covidBVAR:::lag(x, 2)
  expect_equal(z, matrix(c(0, 0, 1, 2, 3), ncol = 1))

  # Test case 3: Custom initial value
  x <- matrix(1:5, ncol = 1)
  z <- covidBVAR:::lag(x, 1, 999)
  expect_equal(z, matrix(c(999, 1, 2, 3, 4), ncol = 1))

  # Test case 4: Matrix with multiple columns
  x <- matrix(1:10, ncol = 2)
  z <- covidBVAR:::lag(x, 1)
  expect_equal(nrow(z), 5)
  expect_equal(ncol(z), 2)
  expect_equal(z[1, ], c(0, 0))

  # Test case 5: n <= 0 returns NULL
  x <- matrix(1:5, ncol = 1)
  z <- covidBVAR:::lag(x, 0)
  expect_null(z)

  z <- covidBVAR:::lag(x, -1)
  expect_null(z)
})


test_that("ols1 works correctly", {
  # Test case: Simple linear regression
  set.seed(123)
  n <- 100
  x <- cbind(1, rnorm(n))  # Intercept + one variable
  beta_true <- c(2, 3)
  y <- x %*% beta_true + rnorm(n, sd = 0.5)

  result <- covidBVAR:::ols1(y, x)

  # Check structure
  expect_equal(result$nobs, 100)
  expect_equal(result$nvar, 2)

  # Check that estimates are close to true values
  expect_equal(length(result$bhatols), 2)
  expect_true(abs(result$bhatols[1] - 2) < 0.5)  # Intercept
  expect_true(abs(result$bhatols[2] - 3) < 0.5)  # Slope

  # Check dimensions
  expect_equal(dim(result$bhatols), c(2, 1))
  expect_equal(dim(result$yhatols), c(100, 1))
  expect_equal(dim(result$resols), c(100, 1))

  # Check that residuals are orthogonal to X
  expect_true(all(abs(t(x) %*% result$resols) < 1e-10))

  # Check R2 is between 0 and 1
  expect_true(result$R2 >= 0 && result$R2 <= 1)

  # Verify formula: sig2hatols = (residuals' * residuals) / (n - k)
  expected_sig2 <- as.numeric(t(result$resols) %*% result$resols) / (n - 2)
  expect_equal(result$sig2hatols, expected_sig2)
})


test_that("gamma_coef works correctly", {
  # Test case 1: Mode = 0.2, SD = 0.4
  result <- covidBVAR:::gamma_coef(0.2, 0.4)

  # Check structure
  expect_true("k" %in% names(result))
  expect_true("theta" %in% names(result))

  # Verify the Gamma distribution has correct mode and variance
  # Mode of Gamma(k, theta) = (k-1)*theta (for k > 1)
  # Variance = k * theta^2

  k <- result$k
  theta <- result$theta

  # Check that k > 1 (required for mode formula)
  expect_true(k > 1)

  # Verify mode formula
  calculated_mode <- (k - 1) * theta
  expect_equal(calculated_mode, 0.2, tolerance = 1e-10)

  # Verify variance formula
  calculated_var <- k * theta^2
  expect_equal(sqrt(calculated_var), 0.4, tolerance = 1e-10)

  # Test case 2: Different values
  result2 <- covidBVAR:::gamma_coef(1, 1)
  expect_true(result2$k > 0)
  expect_true(result2$theta > 0)
})


test_that("utility functions handle edge cases", {
  # Empty matrices
  expect_error(covidBVAR:::trimr(matrix(nrow = 0, ncol = 2), 0, 0))

  # Single element
  x <- matrix(5, ncol = 1)
  expect_equal(covidBVAR:::trimr(x, 0, 0), x)
})
