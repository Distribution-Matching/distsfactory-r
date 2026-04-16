# Gamma distribution matching
# Base R parameterization: shape, rate (or scale = 1/rate)
# mean = shape / rate, var = shape / rate^2

gamma_from_mean_var <- function(mu, var) {
  shape <- mu^2 / var
  rate <- mu / var
  new_dist("gamma", list(shape = shape, rate = rate), dgamma, pgamma, qgamma, rgamma)
}

gamma_from_mean_mode <- function(mu, mode) {
  if (mu <= mode) stop(sprintf("Gamma mean (%g) must be greater than mode (%g)", mu, mode))
  theta <- mu - mode  # scale
  shape <- mu / theta
  rate <- 1 / theta
  new_dist("gamma", list(shape = shape, rate = rate), dgamma, pgamma, qgamma, rgamma)
}

gamma_from_mode_var <- function(mode, var) {
  if (mode <= 0) stop(sprintf("Gamma mode must be positive, got %g", mode))
  f <- function(log_alpha) {
    alpha <- exp(log_alpha) + 1  # ensure alpha > 1
    theta <- sqrt(var / alpha)
    (alpha - 1) * theta - mode
  }
  log_alpha <- find_root_1d(f, x0 = 0)
  alpha <- exp(log_alpha) + 1
  theta <- sqrt(var / alpha)
  new_dist("gamma", list(shape = alpha, rate = 1 / theta), dgamma, pgamma, qgamma, rgamma)
}

gamma_from_mode_quantile <- function(mode, p, q) {
  if (mode <= 0) stop(sprintf("Gamma mode must be positive, got %g", mode))
  f <- function(log_alpha) {
    alpha <- exp(log_alpha) + 1
    theta <- mode / (alpha - 1)
    qgamma(p, shape = alpha, scale = theta) - q
  }
  log_alpha <- find_root_1d(f, x0 = 0)
  alpha <- exp(log_alpha) + 1
  theta <- mode / (alpha - 1)
  new_dist("gamma", list(shape = alpha, rate = 1 / theta), dgamma, pgamma, qgamma, rgamma)
}

gamma_from_mode_iqr <- function(mode, iqr) {
  if (mode <= 0) stop(sprintf("Gamma mode must be positive, got %g", mode))
  f <- function(log_alpha) {
    alpha <- exp(log_alpha) + 1
    theta <- mode / (alpha - 1)
    qgamma(0.75, shape = alpha, scale = theta) - qgamma(0.25, shape = alpha, scale = theta) - iqr
  }
  log_alpha <- find_root_1d(f, x0 = 0)
  alpha <- exp(log_alpha) + 1
  theta <- mode / (alpha - 1)
  new_dist("gamma", list(shape = alpha, rate = 1 / theta), dgamma, pgamma, qgamma, rgamma)
}

gamma_from_two_quantiles <- function(p1, q1, p2, q2) {
  r <- q2 / q1
  f <- function(log_alpha) {
    alpha <- exp(log_alpha)
    z1 <- qgamma(p1, shape = alpha, scale = 1)
    z2 <- qgamma(p2, shape = alpha, scale = 1)
    z2 / z1 - r
  }
  log_alpha <- find_root_1d(f, x0 = 0)
  alpha <- exp(log_alpha)
  theta <- q1 / qgamma(p1, shape = alpha, scale = 1)
  new_dist("gamma", list(shape = alpha, rate = 1 / theta), dgamma, pgamma, qgamma, rgamma)
}

gamma_from_mean_quantile <- function(mu, p, q) {
  # mean = alpha * theta, so theta = mu / alpha
  f <- function(log_alpha) {
    alpha <- exp(log_alpha)
    theta <- mu / alpha
    qgamma(p, shape = alpha, scale = theta) - q
  }
  log_alpha <- find_root_1d(f, x0 = 0)
  alpha <- exp(log_alpha)
  theta <- mu / alpha
  new_dist("gamma", list(shape = alpha, rate = 1 / theta), dgamma, pgamma, qgamma, rgamma)
}

gamma_exists_mean_var <- function(mu, var) {
  mu > 0 && var > 0
}

gamma_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = gamma_from_mean_var(spec$mean, spec$var),
    MeanModeSpec = gamma_from_mean_mode(spec$mean, spec$mode),
    ModeVarSpec = gamma_from_mode_var(spec$mode, spec$var),
    ModeQuantileSpec = gamma_from_mode_quantile(spec$mode, spec$p, spec$q),
    ModeIQRSpec = gamma_from_mode_iqr(spec$mode, spec$iqr),
    TwoQuantileSpec = gamma_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    MeanQuantileSpec = gamma_from_mean_quantile(spec$mean, spec$p, spec$q),
    stop(sprintf("Gamma does not support specification type '%s'", cls))
  )
}
