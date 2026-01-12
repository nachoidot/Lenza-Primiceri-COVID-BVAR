# test-priors.R
# Unit tests for prior specification functions

test_that("set_priors returns correct structure", {
  # Call with default parameters
  result <- covidBVAR:::set_priors()

  # Check main structure elements exist
  expect_true("setpriors" %in% names(result))
  expect_true("MIN" %in% names(result))
  expect_true("MAX" %in% names(result))
  expect_true("priorcoef" %in% names(result))
  expect_true("M" %in% names(result))
  expect_true("N" %in% names(result))
  expect_true("const" %in% names(result))

  # Check setpriors sublist
  expect_true("hyperpriors" %in% names(result$setpriors))
  expect_true("Vc" %in% names(result$setpriors))
  expect_true("hz" %in% names(result$setpriors))
})


test_that("set_priors has correct default values", {
  result <- covidBVAR:::set_priors()

  # Check default values match MATLAB setpriors_covid.m
  expect_equal(result$setpriors$hyperpriors, 1)
  expect_equal(result$setpriors$Vc, 10e6)
  expect_equal(result$setpriors$MNalpha, 0)
  expect_equal(result$setpriors$sur, 1)
  expect_equal(result$setpriors$noc, 1)
  expect_equal(result$setpriors$Fcast, 1)
  expect_equal(result$setpriors$mcmc, 0)
  expect_equal(result$setpriors$Ndraws, 20000)
  expect_equal(result$setpriors$MCMCconst, 1)
})


test_that("set_priors expands hz correctly", {
  # Scalar hz should expand to 1:hz
  result <- covidBVAR:::set_priors(hz = 8)
  expect_equal(result$setpriors$hz, 1:8)

  result <- covidBVAR:::set_priors(hz = 12)
  expect_equal(result$setpriors$hz, 1:12)

  # Vector hz should remain unchanged
  result <- covidBVAR:::set_priors(hz = c(1, 4, 8))
  expect_equal(result$setpriors$hz, c(1, 4, 8))
})


test_that("set_priors has correct hyperparameter bounds", {
  result <- covidBVAR:::set_priors()

  # Check MIN bounds (exact MATLAB values)
  expect_equal(result$MIN$lambda, 0.0001)
  expect_equal(result$MIN$alpha, 0.1)
  expect_equal(result$MIN$theta, 0.0001)
  expect_equal(result$MIN$miu, 0.0001)
  expect_equal(result$MIN$eta, c(1, 1, 1, 0.005))

  # Check MAX bounds (exact MATLAB values)
  expect_equal(result$MAX$lambda, 5)
  expect_equal(result$MAX$miu, 50)
  expect_equal(result$MAX$theta, 50)
  expect_equal(result$MAX$alpha, 5)
  expect_equal(result$MAX$eta, c(500, 500, 500, 0.995))
})


test_that("set_priors computes correct hyperprior coefficients", {
  result <- covidBVAR:::set_priors(hyperpriors = 1)

  # Check priorcoef structure exists
  expect_true("lambda" %in% names(result$priorcoef))
  expect_true("miu" %in% names(result$priorcoef))
  expect_true("theta" %in% names(result$priorcoef))
  expect_true("eta4" %in% names(result$priorcoef))

  # Check Gamma coefficients have k and theta
  expect_true("k" %in% names(result$priorcoef$lambda))
  expect_true("theta" %in% names(result$priorcoef$lambda))

  # Check Beta prior on eta4 has alpha and beta
  expect_true("alpha" %in% names(result$priorcoef$eta4))
  expect_true("beta" %in% names(result$priorcoef$eta4))

  # Verify eta4 Beta parameters (should be approximately 3.04, 1.51 for mode=0.8, sd=0.2)
  expect_equal(result$priorcoef$eta4$alpha, 3.04, tolerance = 0.01)
  expect_equal(result$priorcoef$eta4$beta, 1.51, tolerance = 0.01)
})


test_that("set_priors handles hyperpriors=0 correctly", {
  result <- covidBVAR:::set_priors(hyperpriors = 0)

  # priorcoef should be empty list when hyperpriors=0
  expect_equal(length(result$priorcoef), 0)
})


test_that("set_priors handles MCMC parameters correctly", {
  result <- covidBVAR:::set_priors(mcmc = 1, Ndraws = 10000, Ndrawsdiscard = 2500)

  expect_equal(result$M, 10000)
  expect_equal(result$N, 2500)

  # Default Ndrawsdiscard should be Ndraws/2
  result2 <- covidBVAR:::set_priors(mcmc = 1, Ndraws = 10000)
  expect_equal(result2$N, 5000)
})


test_that("beta_coef_solve returns correct values for default case", {
  # Default case: mode = 0.8, sd = 0.2
  result <- covidBVAR:::beta_coef_solve(c(0.8, 0.2))

  # Should return approximately (3.04, 1.51)
  expect_equal(result[1], 3.04, tolerance = 0.01)
  expect_equal(result[2], 1.51, tolerance = 0.01)

  # Verify the solution: mode = (alpha - 1) / (alpha + beta - 2)
  alpha <- result[1]
  beta <- result[2]
  calculated_mode <- (alpha - 1) / (alpha + beta - 2)
  expect_equal(calculated_mode, 0.8, tolerance = 0.001)

  # Verify variance: var = alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
  calculated_var <- alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
  expect_equal(sqrt(calculated_var), 0.2, tolerance = 0.001)
})
