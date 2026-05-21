# Geometric distribution. Base R: dgeom(x, prob).
# Counts failures before first success: P(X=k) = prob*(1-prob)^k.
# mean = (1-p)/p; var = (1-p)/p^2 = mean*(1+mean).
# So p = 1/(1+mean).

geometric_from_mean <- function(mean) {
  if (mean <= 0) stop("Geometric mean must be positive")
  p <- 1 / (1 + mean)
  new_dist("geometric", list(prob = p), dgeom, pgeom, qgeom, rgeom)
}

geometric_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("Geometric mean and var must be positive")
  expected_var <- mean * (1 + mean)
  if (!isTRUE(all.equal(var, expected_var, tolerance = 1e-9)))
    stop(sprintf("Geometric requires var = mean*(1+mean) = %g; got %g",
                 expected_var, var))
  geometric_from_mean(mean)
}

geometric_from_var <- function(var) {
  # var = mean*(1+mean) -> mean^2 + mean - var = 0 -> mean = (-1 + sqrt(1+4*var))/2
  mean <- (-1 + sqrt(1 + 4 * var)) / 2
  geometric_from_mean(mean)
}

geometric_from_quantile <- function(p, q) {
  if (q < 0 || abs(q - round(q)) > 1e-9)
    stop("Geometric quantile target must be a non-negative integer")
  # P(X <= k) = 1 - (1-prob)^(k+1) = p  ->  prob = 1 - (1-p)^(1/(k+1))
  prob <- 1 - (1 - p)^(1 / (q + 1))
  if (!(0 < prob && prob < 1))
    stop("Geometric single-quantile spec degenerate")
  new_dist("geometric", list(prob = prob), dgeom, pgeom, qgeom, rgeom)
}

geometric_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && mean > 0 && var > 0 &&
    isTRUE(all.equal(var, mean * (1 + mean), tolerance = 1e-9))
}

geometric_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanSpec     = geometric_from_mean(spec$mean),
    MeanVarSpec  = geometric_from_mean_var(spec$mean, spec$var),
    VarSpec      = geometric_from_var(spec$var),
    QuantileSpec = geometric_from_quantile(spec$p, spec$q),
    stop(sprintf("Geometric does not support specification type '%s'", cls))
  )
}
