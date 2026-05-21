# Normal distribution matching.
# Base R parameterization: dnorm(x, mean, sd)
# mean = mu, var = sigma^2
# Mode equals mean (symmetric).

normal_from_mean_var <- function(mu, var) {
  if (var <= 0) stop(sprintf("Normal variance must be positive, got %g", var))
  new_dist("normal", list(mean = mu, sd = sqrt(var)),
           dnorm, pnorm, qnorm, rnorm)
}

normal_from_two_quantiles <- function(p1, q1, p2, q2) {
  z1 <- qnorm(p1); z2 <- qnorm(p2)
  if (z2 == z1) stop("normal: two quantiles must specify distinct probabilities")
  sigma <- (q2 - q1) / (z2 - z1)
  if (sigma <= 0) stop(sprintf("normal: implied sigma must be positive, got %g", sigma))
  mu <- q1 - sigma * z1
  new_dist("normal", list(mean = mu, sd = sigma), dnorm, pnorm, qnorm, rnorm)
}

normal_from_mean_quantile <- function(mu, p, q) {
  z <- qnorm(p)
  if (z == 0) stop("normal: quantile at p=0.5 carries no scale information; pin sigma differently")
  sigma <- (q - mu) / z
  if (sigma <= 0) stop(sprintf("normal: implied sigma must be positive, got %g", sigma))
  new_dist("normal", list(mean = mu, sd = sigma), dnorm, pnorm, qnorm, rnorm)
}

# Symmetric: mode == mean.
normal_from_mode_var <- function(mode, var) {
  normal_from_mean_var(mode, var)
}

normal_from_mean_mode <- function(mu, mode) {
  if (!isTRUE(all.equal(mu, mode)))
    stop(sprintf("normal: mean (%g) and mode (%g) must coincide", mu, mode))
  stop("normal: mean+mode alone underdetermines sigma; provide var/std/cv or a quantile")
}

normal_exists_mean_var <- function(mu, var) {
  is.finite(mu) && is.finite(var) && var > 0
}

normal_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec     = normal_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec = normal_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    MeanQuantileSpec = normal_from_mean_quantile(spec$mean, spec$p, spec$q),
    ModeVarSpec     = normal_from_mode_var(spec$mode, spec$var),
    MeanModeSpec    = normal_from_mean_mode(spec$mean, spec$mode),
    stop(sprintf("Normal does not support specification type '%s'", cls))
  )
}
