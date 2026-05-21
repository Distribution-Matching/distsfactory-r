# Erlang distribution = Gamma with integer shape. Parameters: shape (k, int),
# scale (theta > 0). mean = k*theta, var = k*theta^2.
# We wrap base R's gamma functions and surface (shape = k, scale = theta).

erlang_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("Erlang mean and var must be positive")
  k_real <- mean^2 / var
  k <- round(k_real)
  if (k < 1 || abs(k - k_real) > 1e-6)
    stop(sprintf("Erlang mean+var requires integer shape; got %g", k_real))
  theta <- var / mean
  new_dist("erlang", list(shape = k, scale = theta),
           dgamma, pgamma, qgamma, rgamma)
}

erlang_from_mean <- function(mean) {
  # Implicit k=1 (exponential) — but Julia/Python don't use Erlang for k=1.
  # Disallow to keep API explicit; encourage exponential.
  stop("Erlang from mean alone is ambiguous; specify mean+var (with integer k = mean^2/var)")
}

erlang_exists_mean_var <- function(mean, var) {
  if (!is.finite(mean) || !is.finite(var) || mean <= 0 || var <= 0) return(FALSE)
  k_real <- mean^2 / var
  k_real >= 1 - 1e-9 && abs(k_real - round(k_real)) <= 1e-6
}

erlang_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = erlang_from_mean_var(spec$mean, spec$var),
    stop(sprintf("Erlang does not support specification type '%s'", cls))
  )
}
