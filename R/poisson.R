# Poisson distribution. One parameter: lambda > 0. mean = var = lambda.

poisson_from_mean <- function(mean) {
  if (mean <= 0) stop("Poisson mean must be positive")
  new_dist("poisson", list(lambda = mean), dpois, ppois, qpois, rpois)
}

poisson_from_mean_var <- function(mean, var) {
  if (!isTRUE(all.equal(mean, var, tolerance = 1e-9)))
    stop(sprintf("Poisson requires var == mean; got mean=%g, var=%g", mean, var))
  poisson_from_mean(mean)
}

poisson_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && mean > 0 &&
    isTRUE(all.equal(mean, var, tolerance = 1e-9))
}

poisson_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanSpec    = poisson_from_mean(spec$mean),
    MeanVarSpec = poisson_from_mean_var(spec$mean, spec$var),
    stop(sprintf("Poisson does not support specification type '%s'", cls))
  )
}
