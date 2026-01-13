# optimization.R
# Chris Sims' constrained minimization routines
# Exact ports of csminwel.m, csminit.m, numgrad.m, and bfgsi.m

#' Constrained minimization using BFGS
#'
#' Minimizes a function using BFGS algorithm with line search.
#' Exact port of Chris Sims' csminwel.m function.
#'
#' @param fcn Function to minimize (takes x and ... as arguments)
#' @param x0 Initial parameter vector
#' @param H0 Initial inverse Hessian (must be positive definite)
#' @param grad Either NULL (use numerical gradient) or function that computes gradient
#' @param crit Convergence criterion
#' @param nit Maximum number of iterations
#' @param ... Additional arguments passed to fcn
#'
#' @return List with:
#' \describe{
#'   \item{fh}{Function value at minimum}
#'   \item{xh}{Parameter vector at minimum}
#'   \item{gh}{Gradient at minimum}
#'   \item{H}{Inverse Hessian at minimum}
#'   \item{itct}{Iteration count}
#'   \item{fcount}{Function evaluation count}
#'   \item{retcodeh}{Return code}
#' }
#'
#' @note Direct port from MATLAB csminwel.m by Christopher Sims (273 lines)
#'
#' @noRd
csminwel <- function(fcn, x0, H0, grad = NULL, crit = 1e-16, nit = 1000, ...) {

  # Ensure x0 is a column vector
  if (is.vector(x0)) {
    x0 <- matrix(x0, ncol = 1)
  }
  nx <- max(nrow(x0), ncol(x0))

  Verbose <- 1
  NumGrad <- is.null(grad)
  done <- FALSE
  itct <- 0
  fcount <- 0

  # Evaluate function at initial point
  f0 <- fcn(x0, ...)

  if (f0 > 1e50) {
    cat("Bad initial parameter.\n")
    return(list(fh = f0, xh = x0, gh = NULL, H = H0, itct = 0, fcount = 1, retcodeh = -1))
  }

  # Compute initial gradient
  if (NumGrad) {
    grad_result <- numgrad(fcn, x0, ...)
    g <- grad_result$g
    badg <- grad_result$badg
  } else {
    grad_result <- grad(x0, ...)
    g <- grad_result$g
    badg <- grad_result$badg
  }

  retcode3 <- 101
  x <- x0
  f <- f0
  H <- H0
  cliff <- 0

  ## Main optimization loop
  while (!done) {
    g1 <- NULL
    g2 <- NULL
    g3 <- NULL
    badg1 <- 1
    badg2 <- 1
    badg3 <- 1

    cat('-----------------\n')
    cat('-----------------\n')
    cat(sprintf('f at the beginning of new iteration, %20.10f\n', f))
    cat(sprintf('x = %s\n', paste(sprintf('%15.8g', x), collapse = ' ')))

    itct <- itct + 1

    # Line search
    csminit_result <- csminit(fcn = fcn, x0 = x, f0 = f, g0 = g, bad_grad = badg, H0 = H, ...)
    f1 <- csminit_result$fhat
    x1 <- csminit_result$xhat
    fc <- csminit_result$fcount
    retcode1 <- csminit_result$retcode

    fcount <- fcount + fc

    if (retcode1 != 1) {
      if (retcode1 == 2 || retcode1 == 4) {
        wall1 <- 1
        badg1 <- 1
      } else {
        if (NumGrad) {
          grad_result1 <- numgrad(fcn, x1, ...)
          g1 <- grad_result1$g
          badg1 <- grad_result1$badg
        } else {
          grad_result1 <- grad(x1, ...)
          g1 <- grad_result1$g
          badg1 <- grad_result1$badg
        }
        wall1 <- badg1
      }

      if (wall1 && (length(H) > 1)) {
        # Bad gradient or cliff edge - perturb search direction
        Hcliff <- H + diag(diag(H) * runif(nx))
        cat('Cliff. Perturbing search direction.\n')

        csminit_result2 <- csminit(fcn = fcn, x0 = x, f0 = f, g0 = g, bad_grad = badg, H0 = Hcliff, ...)
        f2 <- csminit_result2$fhat
        x2 <- csminit_result2$xhat
        fc <- csminit_result2$fcount
        retcode2 <- csminit_result2$retcode

        fcount <- fcount + fc

        if (f2 < f) {
          if (retcode2 == 2 || retcode2 == 4) {
            wall2 <- 1
            badg2 <- 1
          } else {
            if (NumGrad) {
              grad_result2 <- numgrad(fcn, x2, ...)
              g2 <- grad_result2$g
              badg2 <- grad_result2$badg
            } else {
              grad_result2 <- grad(x2, ...)
              g2 <- grad_result2$g
              badg2 <- grad_result2$badg
            }
            wall2 <- badg2
          }

          if (wall2) {
            cat('Cliff again. Try traversing\n')
            if (norm(x2 - x1, type = "2") < 1e-13) {
              f3 <- f
              x3 <- x
              badg3 <- 1
              retcode3 <- 101
            } else {
              gcliff <- ((f2 - f1) / (norm(x2 - x1, type = "2")^2)) * (x2 - x1)
              if (ncol(x0) > 1) {
                gcliff <- t(gcliff)
              }

              csminit_result3 <- csminit(fcn = fcn, x0 = x, f0 = f, g0 = gcliff, bad_grad = 0, H0 = diag(nx), ...)
              f3 <- csminit_result3$fhat
              x3 <- csminit_result3$xhat
              fc <- csminit_result3$fcount
              retcode3 <- csminit_result3$retcode

              fcount <- fcount + fc

              if (retcode3 == 2 || retcode3 == 4) {
                wall3 <- 1
                badg3 <- 1
              } else {
                if (NumGrad) {
                  grad_result3 <- numgrad(fcn, x3, ...)
                  g3 <- grad_result3$g
                  badg3 <- grad_result3$badg
                } else {
                  grad_result3 <- grad(x3, ...)
                  g3 <- grad_result3$g
                  badg3 <- grad_result3$badg
                }
                wall3 <- badg3
              }
            }
          } else {
            f3 <- f
            x3 <- x
            badg3 <- 1
            retcode3 <- 101
          }
        } else {
          f3 <- f
          x3 <- x
          badg3 <- 1
          retcode3 <- 101
        }
      } else {
        # Normal iteration
        f2 <- f
        f3 <- f
        badg2 <- 1
        badg3 <- 1
        retcode2 <- 101
        retcode3 <- 101
      }
    } else {
      f2 <- f
      f3 <- f
      f1 <- f
      retcode2 <- retcode1
      retcode3 <- retcode1
    }

    ## Pick best result
    if (f3 < f - crit && badg3 == 0) {
      ih <- 3
      fh <- f3
      xh <- x3
      gh <- g3
      badgh <- badg3
      retcodeh <- retcode3
    } else if (f2 < f - crit && badg2 == 0) {
      ih <- 2
      fh <- f2
      xh <- x2
      gh <- g2
      badgh <- badg2
      retcodeh <- retcode2
    } else if (f1 < f - crit && badg1 == 0) {
      ih <- 1
      fh <- f1
      xh <- x1
      gh <- g1
      badgh <- badg1
      retcodeh <- retcode1
    } else {
      fvals <- c(f1, f2, f3)
      ih <- which.min(fvals)
      fh <- fvals[ih]
      cat(sprintf('ih = %d\n', ih))

      if (ih == 1) {
        xh <- x1
      } else if (ih == 2) {
        xh <- x2
      } else {
        xh <- x3
      }

      retcodei <- c(retcode1, retcode2, retcode3)
      retcodeh <- retcodei[ih]

      # Compute gradient if needed
      if (NumGrad) {
        grad_result_h <- numgrad(fcn, xh, ...)
        gh <- grad_result_h$g
        badgh <- grad_result_h$badg
      } else {
        grad_result_h <- grad(xh, ...)
        gh <- grad_result_h$g
        badgh <- grad_result_h$badg
      }
      badgh <- 1
    }

    ## Check convergence
    stuck <- (abs(fh - f) < crit)

    # Update inverse Hessian (BFGS update)
    if (!badg && !badgh && !stuck) {
      H <- bfgsi(H, gh - g, xh - x)
    }

    if (Verbose) {
      cat('----\n')
      cat(sprintf('Improvement on iteration %d = %18.9f\n', itct, f - fh))
    }

    # Check termination conditions
    if (itct > nit) {
      cat('iteration count termination\n')
      done <- TRUE
    } else if (stuck) {
      cat('improvement < crit termination\n')
      done <- TRUE
    }

    rc <- retcodeh
    if (rc == 1) {
      cat('zero gradient\n')
    } else if (rc == 6) {
      cat('smallest step still improving too slow, reversed gradient\n')
    } else if (rc == 5) {
      cat('largest step still improving too fast\n')
    } else if (rc == 4 || rc == 2) {
      cat('back and forth on step length never finished\n')
    } else if (rc == 3) {
      cat('smallest step still improving too slow\n')
    } else if (rc == 7) {
      cat('warning: possible inaccuracy in H matrix\n')
    }

    # Update for next iteration
    f <- fh
    x <- xh
    g <- gh
    badg <- badgh
  }

  return(list(fh = fh, xh = xh, gh = gh, H = H, itct = itct, fcount = fcount, retcodeh = retcodeh))
}


#' Line search initialization for csminwel
#'
#' Performs line search along descent direction.
#' Exact port of Chris Sims' csminit.m function.
#'
#' @note Direct port from MATLAB csminit.m (193 lines)
#'
#' @noRd
csminit <- function(fcn, x0, f0, g0, bad_grad, H0, ...) {
  # NOTE: Parameter renamed from 'badg' to 'bad_grad' to avoid R's partial
  # argument matching with 'b' (prior mean matrix) passed through ...

  ANGLE <- 0.005
  THETA <- 0.3
  FCHANGE <- 1000
  MINLAMB <- 1e-9
  MINDFAC <- 0.01

  fcount <- 0
  lambda <- 1
  xhat <- x0
  f <- f0
  fhat <- f0
  g <- g0
  gnorm <- norm(g, type = "2")

  if (gnorm < 1e-12 && !bad_grad) {
    retcode <- 1
    dxnorm <- 0
    # Gradient convergence
  } else {
    # Gauss-Newton step
    dx <- -H0 %*% g
    dxnorm <- norm(dx, type = "2")

    if (dxnorm > 1e12) {
      cat('Near-singular H problem.\n')
      dx <- dx * FCHANGE / dxnorm
    }

    dfhat <- drop(t(dx) %*% g0)

    if (!bad_grad) {
      # Test for alignment of dx with gradient
      a <- -dfhat / (gnorm * dxnorm)
      if (a < ANGLE) {
        dx <- dx - (ANGLE * dxnorm / gnorm + dfhat / (gnorm * gnorm)) * g
        dx <- dx * dxnorm / norm(dx, type = "2")  # Keep scale invariant
        dfhat <- drop(t(dx) %*% g)
        cat(sprintf('Correct for low angle: %g\n', a))
      }
    }

    cat(sprintf('Predicted improvement: %18.9f\n', -dfhat / 2))

    # Adjust step length (lambda)
    done <- FALSE
    factor <- 3
    shrink <- 1
    lambdaMin <- 0
    lambdaMax <- Inf
    lambdaPeak <- 0
    fPeak <- f0
    lambdahat <- 0

    while (!done) {
      if (ncol(x0) > 1) {
        dxtest <- x0 + t(dx) * lambda
      } else {
        dxtest <- x0 + dx * lambda
      }

      f <- fcn(dxtest, ...)
      cat(sprintf('lambda = %10.5g; f = %20.7f\n', lambda, f))

      # Handle NaN/Inf like MATLAB: treat as bad step, try smaller lambda
      if (is.nan(f) || is.infinite(f)) {
        shrinkSignal <- TRUE
        growSignal <- FALSE
      } else {
        if (f < fhat) {
          fhat <- f
          xhat <- dxtest
          lambdahat <- lambda
        }

        shrinkSignal <- (!bad_grad && (f0 - f < max(c(-THETA * dfhat * lambda, 0)))) ||
                        (bad_grad && (f0 - f) < 0)
        growSignal <- !bad_grad && (lambda > 0) && (f0 - f > -(1 - THETA) * dfhat * lambda)
      }

      fcount <- fcount + 1

      if (shrinkSignal && ((lambda > lambdaPeak) || (lambda < 0))) {
        if ((lambda > 0) && ((!shrink) || (lambda / factor <= lambdaPeak))) {
          shrink <- 1
          factor <- factor^0.6
          while (lambda / factor <= lambdaPeak) {
            factor <- factor^0.6
          }
          if (abs(factor - 1) < MINDFAC) {
            if (abs(lambda) < 4) {
              retcode <- 2
            } else {
              retcode <- 7
            }
            done <- TRUE
          }
        }
        if ((lambda < lambdaMax) && (lambda > lambdaPeak)) {
          lambdaMax <- lambda
        }
        lambda <- lambda / factor
        if (abs(lambda) < MINLAMB) {
          if ((lambda > 0) && (f0 <= fhat)) {
            lambda <- -lambda * factor^6
          } else {
            if (lambda < 0) {
              retcode <- 6
            } else {
              retcode <- 3
            }
            done <- TRUE
          }
        }
      } else if ((growSignal && lambda > 0) ||
                 (shrinkSignal && ((lambda <= lambdaPeak) && (lambda > 0)))) {
        if (shrink) {
          shrink <- 0
          factor <- factor^0.6
          if (abs(factor - 1) < MINDFAC) {
            if (abs(lambda) < 4) {
              retcode <- 4
            } else {
              retcode <- 7
            }
            done <- TRUE
          }
        }
        if ((f < fPeak) && (lambda > 0)) {
          fPeak <- f
          lambdaPeak <- lambda
          if (lambdaMax <= lambdaPeak) {
            lambdaMax <- lambdaPeak * factor * factor
          }
        }
        lambda <- lambda * factor
        if (abs(lambda) > 1e20) {
          retcode <- 5
          done <- TRUE
        }
      } else {
        done <- TRUE
        if (factor < 1.2) {
          retcode <- 7
        } else {
          retcode <- 0
        }
      }
    }
  }

  cat(sprintf('Norm of dx %10.5g\n', dxnorm))

  return(list(fhat = fhat, xhat = xhat, fcount = fcount, retcode = retcode))
}


#' Numerical gradient
#'
#' Computes numerical gradient using forward differences.
#' Exact port of Chris Sims' numgrad.m function.
#'
#' @param fcn Function to differentiate
#' @param x Point at which to compute gradient
#' @param ... Additional arguments to fcn
#'
#' @return List with:
#' \describe{
#'   \item{g}{Gradient vector}
#'   \item{badg}{Flag indicating bad gradient (0 = good, 1 = bad)}
#' }
#'
#' @note Direct port from MATLAB numgrad.m (97 lines)
#'
#' @noRd
numgrad <- function(fcn, x0_ng, ...) {
  # NOTE: Parameter renamed from 'x' to 'x0_ng' to avoid collision with
  # 'x' (data matrix) that may be passed through ... from bvar_estimate
  delta <- 1e-6
  n <- length(x0_ng)
  tvec <- delta * diag(n)
  g <- rep(0, n)

  f0 <- fcn(x0_ng, ...)
  badg <- 0

  for (i in 1:n) {
    if (nrow(x0_ng) > ncol(x0_ng)) {
      tvecv <- tvec[i, ]
    } else {
      tvecv <- tvec[, i]
    }

    g0 <- (fcn(x0_ng + tvecv, ...) - f0) / delta

    if (abs(g0) < 1e15) {
      g[i] <- g0
    } else {
      cat('bad gradient ------------------------\n')
      g[i] <- 0
      badg <- 1
    }
  }

  return(list(g = g, badg = badg))
}


#' BFGS inverse Hessian update
#'
#' Updates inverse Hessian using BFGS formula.
#' Exact port of Chris Sims' bfgsi.m function.
#'
#' @param H0 Current inverse Hessian
#' @param dg Change in gradient
#' @param dx Change in x
#'
#' @return Updated inverse Hessian
#'
#' @note Direct port from MATLAB bfgsi.m (25 lines)
#'
#' @noRd
bfgsi <- function(H0, dg, dx) {
  # Ensure column vectors
  # Handle both vector and matrix inputs
  if (is.matrix(dg) && ncol(dg) > 1) {
    dg <- t(dg)
  }
  if (!is.matrix(dg)) {
    dg <- matrix(dg, ncol = 1)
  }
  if (is.matrix(dx) && ncol(dx) > 1) {
    dx <- t(dx)
  }
  if (!is.matrix(dx)) {
    dx <- matrix(dx, ncol = 1)
  }

  Hdg <- H0 %*% dg
  dgdx <- drop(t(dg) %*% dx)  # Convert to scalar

  if (abs(dgdx) > 1e-12) {
    dg_Hdg <- drop(t(dg) %*% Hdg)  # Convert to scalar
    H <- H0 + (1 + dg_Hdg / dgdx) * (dx %*% t(dx)) / dgdx -
         (dx %*% t(Hdg) + Hdg %*% t(dx)) / dgdx
  } else {
    cat('bfgs update failed.\n')
    cat(sprintf('|dg| = %g |dx| = %g\n', sqrt(drop(t(dg) %*% dg)), sqrt(drop(t(dx) %*% dx))))
    cat(sprintf("dg'*dx = %g\n", dgdx))
    cat(sprintf('|H*dg| = %g\n', drop(t(Hdg) %*% Hdg)))
    H <- H0
  }

  return(H)
}
