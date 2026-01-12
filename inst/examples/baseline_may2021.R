# baseline_may2021.R
# Replication of Baseline_May2021.m from Lenza & Primiceri (2022)
# Estimates BVAR with COVID volatility using data through May 2021

library(covidBVAR)

## Load data
data(macro_data_may2021)

# Extract baseline model variables (matching MATLAB line 32)
# indmacro=[15 8 9 14 12 6 13]
y_lev <- macro_data_may2021$data[, macro_data_may2021$baseline_indices]
colnames(y_lev) <- macro_data_may2021$baseline_names

# Apply log transformation (MATLAB line 56)
y_log <- 100 * log(y_lev)

# Extract key dates
T0 <- macro_data_may2021$key_dates$T0          # Dec 1988
T1estim <- macro_data_may2021$key_dates$T1estim  # May 2021
Tcovid <- macro_data_may2021$key_dates$Tcovid   # March 2020 (relative to T0)

# Subset to estimation sample
y <- y_log[T0:T1estim, ]
n <- ncol(y)

cat("Estimation sample:", nrow(y), "observations\n")
cat("Variables:", n, "\n")
cat("COVID break at observation:", Tcovid, "\n\n")

## Estimate BVAR (matching MATLAB line 63)
# rng(10); lags=13; ndraws=2*2500;
# res = bvarGLP_covid(Ylog,lags,'mcmc',1,'MCMCconst',1,
#                     'MNpsi',0,'sur',0,'noc',0,'Ndraws',ndraws,
#                     'hyperpriors',1,'Tcovid',Tcovid);

set.seed(10)           # Match MATLAB rng(10)
lags <- 13
ndraws <- 2 * 2500

cat("Starting BVAR estimation with COVID volatility...\n")
cat("----------------------------------------------\n")
cat("Lags:", lags, "\n")
cat("MCMC draws:", ndraws, "\n")
cat("Hyperpriors: Yes (reference priors)\n")
cat("SUR prior: No\n")
cat("NOC prior: No\n")
cat("MNalpha: Fixed at 2\n\n")

res <- bvar_estimate(
  y = y,
  lags = lags,
  hyperpriors = 1,     # Reference priors on hyperparameters
  Tcovid = Tcovid,     # March 2020 volatility break
  sur = 0,             # No sum-of-coefficients prior
  noc = 0,             # No no-cointegration prior
  MNalpha = 0,         # Fix lag decay at 2
  mcmc = 1,            # Run MCMC
  Ndraws = ndraws,
  MCMCconst = 1,       # MCMC proposal scaling
  MCMCfcast = 1,       # Generate forecasts in MCMC
  MCMCstorecoeff = 1   # Store VAR coefficients
)

cat("\n\nEstimation complete!\n")
cat("-------------------\n")

## Display posterior mode estimates
cat("\nPosterior mode hyperparameters:\n")
cat("  Lambda (Minnesota shrinkage):", round(res$postmax$lambda, 4), "\n")
cat("  COVID volatility parameters:\n")
cat("    eta1 (March 2020):", round(res$postmax$eta[1], 2), "\n")
cat("    eta2 (April 2020):", round(res$postmax$eta[2], 2), "\n")
cat("    eta3 (asymptotic):", round(res$postmax$eta[3], 2), "\n")
cat("    eta4 (decay rate):", round(res$postmax$eta[4], 4), "\n")
cat("  Log-posterior:", round(res$postmax$logPost, 2), "\n")
cat("  Optimization iterations:", res$postmax$itct, "\n")

## Display MCMC diagnostics
if (!is.null(res$mcmc)) {
  cat("\nMCMC diagnostics:\n")
  cat("  Acceptance rate:", round(res$mcmc$ACCrate * 100, 1), "%\n")
  cat("  Target: ~25%\n")

  cat("\nPosterior means (MCMC):\n")
  cat("  Lambda:", round(mean(res$mcmc$lambda), 4),
      "(SD:", round(sd(res$mcmc$lambda), 4), ")\n")
  cat("  eta1:", round(mean(res$mcmc$eta[, 1]), 2),
      "(SD:", round(sd(res$mcmc$eta[, 1]), 2), ")\n")
  cat("  eta2:", round(mean(res$mcmc$eta[, 2]), 2),
      "(SD:", round(sd(res$mcmc$eta[, 2]), 2), ")\n")
  cat("  eta3:", round(mean(res$mcmc$eta[, 3]), 2),
      "(SD:", round(sd(res$mcmc$eta[, 3]), 2), ")\n")
  cat("  eta4:", round(mean(res$mcmc$eta[, 4]), 4),
      "(SD:", round(sd(res$mcmc$eta[, 4]), 4), ")\n")
}

## Compute IRFs (matching MATLAB lines 66-73)
# Impulse responses to unemployment shock
H <- 60  # Horizon
if (!is.null(res$mcmc)) {
  M <- dim(res$mcmc$beta)[3]
  Dirf1 <- array(0, dim = c(H + 1, n, M))

  cat("\nComputing impulse response functions...\n")
  for (jg in 1:M) {
    Dirf1[, , jg] <- bvar_irfs(
      beta = res$mcmc$beta[, , jg],
      sigma = res$mcmc$sigma[, , jg],
      nshock = 1,  # Unemployment shock
      hmax = H + 1
    )
  }

  # Sort for percentiles
  sIRF1 <- apply(Dirf1, c(1, 2), sort)

  cat("IRFs computed for", M, "MCMC draws\n")
  cat("Horizon:", H + 1, "periods\n")
}

## Save results
cat("\nSaving results...\n")
save(res, file = "results/Baseline_May2021.RData")
if (exists("Dirf1")) {
  save(Dirf1, sIRF1, file = "results/Baseline_May2021_IRFs.RData")
}

cat("\nDone! Results saved to results/\n")
cat("- Baseline_May2021.RData: Full estimation results\n")
if (exists("Dirf1")) {
  cat("- Baseline_May2021_IRFs.RData: Impulse response functions\n")
}

## Comparison with MATLAB
cat("\n\nNOTE: Compare these estimates with MATLAB output:\n")
cat("------------------------------------------------------\n")
cat("MATLAB results should be in: Example/lp-codes/results/Baseline_May2021.mat\n")
cat("Load in R using: R.matlab::readMat('Baseline_May2021.mat')\n")
cat("Key values to compare:\n")
cat("  - res$postmax$lambda\n")
cat("  - res$postmax$eta\n")
cat("  - res$postmax$logPost\n")
cat("  - res$mcmc$ACCrate\n")
cat("Expected tolerance: < 1e-6 for posterior mode\n")
