# Uniform distribution matching.
# Base R parameterization: dunif(x, min, max)
# mean = (a+b)/2, var = (b-a)^2 / 12
# Mode is not unique (constant density), so mode-based specs are unsupported.

uniform_from_mean_var <- function(mean, var) {
  if (var <= 0) stop(sprintf("Uniform variance must be positive, got %g", var))
  half <- sqrt(3 * var)
  new_dist("uniform", list(min = mean - half, max = mean + half),
           dunif, punif, qunif, runif)
}

uniform_from_two_quantiles <- function(p1, q1, p2, q2) {
  if (p2 == p1) stop("uniform: two-quantile probabilities must differ")
  width <- (q2 - q1) / (p2 - p1)
  if (width <= 0) stop("uniform: implied width must be positive")
  a <- q1 - p1 * width
  b <- a + width
  new_dist("uniform", list(min = a, max = b), dunif, punif, qunif, runif)
}

uniform_from_mean_quantile <- function(mean, p, q) {
  # mean = (a+b)/2; q = a + p*(b-a)
  # Let w = b-a > 0, then a = mean - w/2, b = mean + w/2, and
  # q = (mean - w/2) + p*w = mean + w*(p - 1/2). Solve w = (q - mean)/(p - 1/2).
  if (p == 0.5) stop("uniform: mean+median underdetermines width")
  w <- (q - mean) / (p - 0.5)
  if (w <= 0) stop("uniform: implied width must be positive")
  new_dist("uniform", list(min = mean - w / 2, max = mean + w / 2),
           dunif, punif, qunif, runif)
}

uniform_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && var > 0
}

uniform_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec      = uniform_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec  = uniform_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    MeanQuantileSpec = uniform_from_mean_quantile(spec$mean, spec$p, spec$q),
    stop(sprintf("Uniform does not support specification type '%s'", cls))
  )
}
