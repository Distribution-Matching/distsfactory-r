# Numerical solvers for parameter fitting

#' Find root of a scalar function with auto-bracketing
#' @param f Function of one variable.
#' @param x0 Starting point for bracket search (default 0).
#' @param bracket Optional length-2 vector \code{c(a, b)} where \code{f(a)}
#'   and \code{f(b)} have opposite signs.
#' @return The root.
#' @keywords internal
find_root_1d <- function(f, x0 = 0, bracket = NULL) {
  if (!is.null(bracket)) {
    return(uniroot(f, bracket, tol = 1e-12)$root)
  }

  # Auto-bracket: expand geometrically from x0
  a <- x0 - 1
  b <- x0 + 1
  for (i in seq_len(60)) {
    fa <- tryCatch(f(a), error = function(e) NA_real_)
    fb <- tryCatch(f(b), error = function(e) NA_real_)
    if (is.finite(fa) && is.finite(fb) && fa * fb < 0) {
      return(uniroot(f, c(a, b), tol = 1e-12)$root)
    }
    a <- a * 2
    b <- b * 2
  }
  stop(sprintf("Could not bracket a root starting from x0=%g", x0))
}


#' Damped Newton iteration for a 2D system F(x) = 0
#' @param F Function taking length-2 vector, returning length-2 vector.
#' @param x0 Initial guess (length-2 vector).
#' @param maxiter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param h Finite-difference step size.
#' @return Solution vector (length 2).
#' @keywords internal
newton_2d <- function(F, x0, maxiter = 200, tol = 1e-10, h = 1e-7) {
  x <- x0
  for (iter in seq_len(maxiter)) {
    Fx <- F(x)
    if (max(abs(Fx)) < tol) return(x)

    # Numerical Jacobian
    J <- matrix(0, 2, 2)
    for (j in 1:2) {
      xp <- x
      xp[j] <- xp[j] + h
      J[, j] <- (F(xp) - Fx) / h
    }

    dx <- tryCatch(solve(J, Fx), error = function(e) {
      stop("Singular Jacobian in Newton iteration")
    })

    # Damped step
    step <- 1.0
    for (k in seq_len(20)) {
      x_new <- x - step * dx
      Fx_new <- F(x_new)
      if (max(abs(Fx_new)) < max(abs(Fx))) break
      step <- step * 0.5
    }
    x <- x - step * dx
  }
  stop(sprintf(
    "Newton iteration did not converge after %d iterations (residual: %g)",
    maxiter, max(abs(F(x)))
  ))
}
