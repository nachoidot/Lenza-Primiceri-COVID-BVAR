# bvar_estimate.R
# Main BVAR estimation function
# Exact port of bvarGLP_covid.m (413 lines)

#' Estimate Bayesian VAR with COVID-19 volatility
#'
#' Estimates the BVAR of Giannone, Lenza and Primiceri (2015), augmented with
#' a change in volatility at the time of COVID (March 2020). This is the main
#' user-facing function that integrates all components.
#'
#' @param y Data matrix (T x n) where T is number of observations, n is number of variables
#' @param lags Number of lags in the VAR
#' @param hyperpriors 0 = no priors on hyperparameters, 1 = reference priors (default: 1)
#' @param Vc Prior variance on constant term (default: 10e6)
#' @param pos Variables entering in first differences (default: integer(0))
#' @param MNalpha 0 = fix lag decay at 2, 1 = estimate (default: 0)
#' @param Tcovid Position of March 2020 observation, or NULL for no COVID break (default: NULL)
#' @param sur 0 = SUR prior OFF, 1 = SUR prior ON (default: 1)
#' @param noc 0 = NOC prior OFF, 1 = NOC prior ON (default: 1)
#' @param Fcast 0 = no forecasts, 1 = generate forecasts at posterior mode (default: 1)
#' @param hz Forecast horizons (default: 1:8)
#' @param mcmc 0 = no MCMC, 1 = run MCMC (default: 0)
#' @param Ndraws Number of MCMC draws (default: 20000)
#' @param Ndrawsdiscard Burn-in draws (default: Ndraws/2)
#' @param MCMCconst MCMC proposal scaling (default: 1)
#' @param MCMCfcast Generate forecasts in MCMC (default: 1)
#' @param MCMCstorecoeff Store VAR coefficients (default: 1)
#'
#' @return List with estimation results:
#' \describe{
#'   \item{lags}{Number of lags}
#'   \item{Tcovid}{COVID break position}
#'   \item{postmax}{Posterior mode estimates (betahat, sigmahat, hyperparameters, logPost, etc.)}
#'   \item{mcmc}{MCMC results (if mcmc=1): draws of hyperparameters, coefficients, forecasts}
#'   \item{setpriors}{Prior settings used}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(macro_data_may2021)
#' y <- macro_data_may2021$data
#'
#' # Estimate BVAR with COVID volatility
#' result <- bvar_estimate(y, lags = 13, Tcovid = 376, mcmc = 1, Ndraws = 5000)
#'
#' # Extract posterior mode
#' result$postmax$lambda  # Minnesota shrinkage
#' result$postmax$eta     # COVID volatility parameters
#' }
#'
#' @references
#' Lenza, M., & Primiceri, G. E. (2022). How to estimate a VAR after March 2020.
#' Journal of Applied Econometrics, 37(4), 688-699.
#'
#' @note Direct port from MATLAB bvarGLP_covid.m (413 lines)
bvar_estimate <- function(y,
                          lags,
                          hyperpriors = 1,
                          Vc = 10e6,
                          pos = integer(0),
                          MNalpha = 0,
                          Tcovid = NULL,
                          sur = 1,
                          noc = 1,
                          Fcast = 1,
                          hz = 8,
                          mcmc = 0,
                          Ndraws = 20000,
                          Ndrawsdiscard = NULL,
                          MCMCconst = 1,
                          MCMCfcast = 1,
                          MCMCstorecoeff = 1) {

  ## Set BVAR priors
  prior_settings <- set_priors(
    hyperpriors = hyperpriors,
    Vc = Vc,
    pos = pos,
    MNalpha = MNalpha,
    Tcovid = Tcovid,
    sur = sur,
    noc = noc,
    Fcast = Fcast,
    hz = hz,
    mcmc = mcmc,
    Ndraws = Ndraws,
    Ndrawsdiscard = Ndrawsdiscard,
    MCMCconst = MCMCconst,
    MCMCfcast = MCMCfcast,
    MCMCstorecoeff = MCMCstorecoeff
  )

  MIN <- prior_settings$MIN
  MAX <- prior_settings$MAX
  priorcoef <- prior_settings$priorcoef
  M <- prior_settings$M
  N <- prior_settings$N
  const <- prior_settings$const

  # Use expanded hz from set_priors (converts scalar 8 to 1:8)
  hz <- prior_settings$setpriors$hz

  ## Data matrix manipulations (lines 24-43 in MATLAB)
  if (is.vector(y)) {
    y <- matrix(y, ncol = 1)
  }

  TT <- nrow(y)
  n <- ncol(y)
  k <- n * lags + 1  # Number of coefficients per equation

  # Construct matrix of regressors
  xdata <- matrix(0, nrow = TT, ncol = k)
  xdata[, 1] <- 1  # Constant

  for (i in 1:lags) {
    xdata[, (1 + (i - 1) * n + 1):(1 + i * n)] <- lag(y, i)
  }

  # Initial condition (mean of first lags observations)
  y0 <- colMeans(y[1:lags, , drop = FALSE])

  # Trim initial lags
  xdata <- xdata[(lags + 1):TT, , drop = FALSE]
  y <- y[(lags + 1):TT, , drop = FALSE]
  T <- nrow(y)

  # Adjust Tcovid
  if (!is.null(Tcovid)) {
    Tcovid <- Tcovid - lags
  }

  ## Minnesota prior mean (lines 46-51)
  b <- matrix(0, nrow = k, ncol = n)
  diagb <- rep(1, n)
  diagb[pos] <- 0  # Zero for variables in first differences
  b[2:(n + 1), ] <- diag(diagb)

  ## Starting values for minimization (lines 54-68)
  lambda0 <- 0.2
  theta0 <- 1
  miu0 <- 1
  alpha0 <- 2

  # COVID volatility starting values
  # MATLAB: aux=mean(abs(y(Tcovid:...)' ,1) / mean(mean(abs(y(2:Tcovid-1,:)-y(1:Tcovid-2,:))))
  # The transpose and mean(...,1) in MATLAB gives row-wise means in R terms
  # aux should have length = number of COVID observations, not number of variables
  if (!is.null(Tcovid)) {
    Tend <- min(max(Tcovid + 1, T), T)
    # Row means (mean across variables for each time point)
    # This matches MATLAB's mean(X', 1) which computes mean across variables
    aux <- rowMeans(abs(y[Tcovid:Tend, , drop = FALSE] -
                        y[(Tcovid - 1):(Tend - 1), , drop = FALSE])) /
           mean(abs(y[2:(Tcovid - 1), ] - y[1:(Tcovid - 2), ]))

    if (length(aux) == 0) {
      eta0 <- numeric(0)
    } else if (length(aux) == 2) {
      eta0 <- c(aux, aux[1], 0.8)
    } else if (length(aux) >= 3) {
      eta0 <- c(aux[1:3], 0.8)
    } else {
      eta0 <- c(aux[1], aux[1], aux[1], 0.8)
    }
  } else {
    eta0 <- numeric(0)
  }

  ## Residual variance of AR(1) for each variable (lines 71-77)
  SS <- numeric(n)
  for (i in 1:n) {
    Tend <- T
    if (!is.null(Tcovid)) {
      Tend <- min(Tcovid - 1, T)
    }
    if (Tend > 1) {
      ar1 <- ols1(y[2:Tend, i, drop = FALSE],
                  cbind(1, y[1:(Tend - 1), i, drop = FALSE]))
      SS[i] <- ar1$sig2hatols
    } else {
      SS[i] <- var(y[, i])
    }
  }

  ## Logit transformation for starting values (lines 79-116)
  inlambda <- -log((MAX$lambda - lambda0) / (lambda0 - MIN$lambda))
  inHlambda <- (1 / (MAX$lambda - lambda0) + 1 / (lambda0 - MIN$lambda))^2 *
               (abs(lambda0) / 1)^2

  if (MNalpha == 1) {
    inalpha <- -log((MAX$alpha - alpha0) / (alpha0 - MIN$alpha))
    inHalpha <- (1 / (MAX$alpha - alpha0) + 1 / (alpha0 - MIN$alpha))^2 *
                (abs(alpha0) / 1)^2
  } else {
    inalpha <- numeric(0)
    inHalpha <- numeric(0)
  }

  if (sur == 1) {
    intheta <- -log((MAX$theta - theta0) / (theta0 - MIN$theta))
    inHtheta <- (1 / (MAX$theta - theta0) + 1 / (theta0 - MIN$theta))^2 *
                (abs(theta0) / 1)^2
  } else {
    intheta <- numeric(0)
    inHtheta <- numeric(0)
  }

  if (noc == 1) {
    inmiu <- -log((MAX$miu - miu0) / (miu0 - MIN$miu))
    inHmiu <- (1 / (MAX$miu - miu0) + 1 / (miu0 - MIN$miu))^2 *
              (abs(miu0) / 1)^2
  } else {
    inmiu <- numeric(0)
    inHmiu <- numeric(0)
  }

  if (!is.null(Tcovid) && length(eta0) > 0) {
    ncp <- length(eta0)
    ineta <- -log((MAX$eta - eta0) / (eta0 - MIN$eta))
    inHeta <- (1 / (MAX$eta - eta0) + 1 / (eta0 - MIN$eta))^2 *
              (abs(eta0) / 1)^2
  } else {
    ineta <- numeric(0)
    inHeta <- numeric(0)
  }

  x0 <- c(inlambda, ineta, intheta, inmiu, inalpha)
  H0 <- diag(c(inHlambda, inHeta, inHtheta, inHmiu, inHalpha))

  ## Maximization of posterior (line 120)
  cat('Starting hyperparameter optimization...\n')

  mn <- list(alpha = MNalpha)

  opt_result <- csminwel(
    fcn = function(par, ...) {
      result <- tryCatch(
        log_ml_var_formin(par, ...),
        error = function(e) list(logML = NaN)  # Return NaN on error (matches MATLAB behavior)
      )
      return(result$logML)
    },
    x0 = x0,
    H0 = H0,
    grad = NULL,
    crit = 1e-16,
    nit = 1000,
    y = y,
    xdata = xdata,
    lags = lags,
    T = T,
    n = n,
    b = b,
    MIN = MIN,
    MAX = MAX,
    SS = SS,
    Vc = Vc,
    pos = pos,
    mn = mn,
    sur = sur,
    noc = noc,
    y0 = y0,
    hyperpriors = hyperpriors,
    priorcoef = priorcoef,
    Tcovid = Tcovid
  )

  fh <- opt_result$fh
  xh <- opt_result$xh
  gh <- opt_result$gh
  H <- opt_result$H
  itct <- opt_result$itct
  fcount <- opt_result$fcount
  retcodeh <- opt_result$retcodeh

  ## Output of maximization (lines 124-178)
  postmax_result <- log_ml_var_formin(
    xh, y, xdata, lags, T, n, b, MIN, MAX, SS, Vc, pos,
    mn, sur, noc, y0, hyperpriors, priorcoef, Tcovid
  )

  betahat <- postmax_result$betahat
  sigmahat <- postmax_result$sigmahat

  # Extract hyperparameters from transformed scale
  lambda <- MIN$lambda + (MAX$lambda - MIN$lambda) / (1 + exp(-xh[1]))

  theta <- MAX$theta
  miu <- MAX$miu
  eta <- MAX$eta
  alpha <- 2

  if (!is.null(Tcovid) && length(ineta) > 0) {
    ncp <- length(ineta)
    eta <- MIN$eta + (MAX$eta - MIN$eta) / (1 + exp(-xh[2:(ncp + 1)]))

    if (sur == 1) {
      theta <- MIN$theta + (MAX$theta - MIN$theta) / (1 + exp(-xh[ncp + 2]))
      if (noc == 1) {
        miu <- MIN$miu + (MAX$miu - MIN$miu) / (1 + exp(-xh[ncp + 3]))
      }
    } else if (sur == 0) {
      if (noc == 1) {
        miu <- MIN$miu + (MAX$miu - MIN$miu) / (1 + exp(-xh[ncp + 2]))
      }
    }
  } else {
    eta <- c(1, 1, 1, 1)
    if (sur == 1) {
      theta <- MIN$theta + (MAX$theta - MIN$theta) / (1 + exp(-xh[2]))
      if (noc == 1) {
        miu <- MIN$miu + (MAX$miu - MIN$miu) / (1 + exp(-xh[3]))
      }
    } else if (sur == 0) {
      if (noc == 1) {
        miu <- MIN$miu + (MAX$miu - MIN$miu) / (1 + exp(-xh[2]))
      }
    }
  }

  if (MNalpha == 1) {
    alpha <- MIN$alpha + (MAX$alpha - MIN$alpha) / (1 + exp(-xh[length(xh)]))
  }

  ## Forecasts at posterior mode (lines 181-190)
  forecast_mode <- NULL
  if (Fcast == 1) {
    Y <- rbind(y, matrix(0, nrow = max(hz), ncol = n))
    for (tau in 1:max(hz)) {
      xT <- c(1, as.vector(t(Y[(T + tau - 1):(T + tau - lags), , drop = FALSE])))
      xT <- matrix(xT, nrow = 1)
      Y[T + tau, ] <- xT %*% betahat
    }
    forecast_mode <- Y[T + hz, , drop = FALSE]
  }

  ## Posterior mode output
  postmax <- list(
    betahat = betahat,
    sigmahat = sigmahat,
    lambda = lambda,
    theta = theta,
    miu = miu,
    eta = eta,
    alpha = alpha,
    logPost = -fh,
    itct = itct,
    SSar1 = SS,
    forecast = forecast_mode
  )

  ## MCMC (if requested)
  mcmc_result <- NULL
  if (mcmc == 1) {
    cat('\nStarting MCMC sampling...\n')

    mcmc_result <- run_mcmc(
      y = y,
      xdata = xdata,
      lags = lags,
      T = T,
      n = n,
      b = b,
      MIN = MIN,
      MAX = MAX,
      SS = SS,
      Vc = Vc,
      pos = pos,
      mn = mn,
      sur = sur,
      noc = noc,
      y0 = y0,
      hyperpriors = hyperpriors,
      priorcoef = priorcoef,
      Tcovid = Tcovid,
      postmode = postmax,
      xh = xh,
      M = M,
      N = N,
      const = const,
      MCMCfcast = MCMCfcast,
      MCMCstorecoeff = MCMCstorecoeff,
      hz = hz
    )

    # Store HH in postmax (matching MATLAB r.postmax.HH = HH, line 281)
    postmax$HH <- mcmc_result$HH
  }

  ## Return structure matching MATLAB
  result <- list(
    lags = lags,
    Tcovid = Tcovid,
    postmax = postmax,
    mcmc = mcmc_result,
    setpriors = prior_settings$setpriors
  )

  return(result)
}
