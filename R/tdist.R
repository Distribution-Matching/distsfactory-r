# Student t distribution matching.
# Base R parameterization: dt(x, df) with location 0, scale 1.
# E[X] = 0 for df > 1; Var[X] = df / (df - 2) for df > 2.
# Only one free parameter (df), so only certain spec styles are supported.

tdist_from_mean_var <- function(mean, var) {
  if (!isTRUE(all.equal(mean, 0))) {
    stop(sprintf("TDist has mean 0; got mean=%g", mean))
  }
  if (var <= 1) stop(sprintf("TDist variance must exceed 1; got %g", var))
  # var = df / (df - 2)  ->  df = 2 * var / (var - 1)
  df <- 2 * var / (var - 1)
  new_dist("tdist", list(df = df), dt, pt, qt, rt)
}

tdist_from_var <- function(var) tdist_from_mean_var(0, var)

tdist_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && isTRUE(all.equal(mean, 0)) && var > 1
}

tdist_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = tdist_from_mean_var(spec$mean, spec$var),
    VarSpec     = tdist_from_var(spec$var),
    stop(sprintf("TDist does not support specification type '%s'", cls))
  )
}
