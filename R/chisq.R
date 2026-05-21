# Chi-squared distribution matching.
# Base R parameterization: dchisq(x, df)
# mean = df, var = 2*df. Only one free parameter.

chisq_from_mean <- function(mean) {
  if (mean <= 0) stop(sprintf("Chi-squared mean (=df) must be positive, got %g", mean))
  new_dist("chisq", list(df = mean), dchisq, pchisq, qchisq, rchisq)
}

chisq_from_mean_var <- function(mean, var) {
  if (mean <= 0) stop("Chi-squared mean must be positive")
  if (!isTRUE(all.equal(var, 2 * mean))) {
    stop(sprintf("Chi-squared requires var = 2*mean; got mean=%g, var=%g (expected %g)",
                 mean, var, 2 * mean))
  }
  new_dist("chisq", list(df = mean), dchisq, pchisq, qchisq, rchisq)
}

chisq_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && mean > 0 && isTRUE(all.equal(var, 2 * mean))
}

chisq_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanSpec    = chisq_from_mean(spec$mean),
    MeanVarSpec = chisq_from_mean_var(spec$mean, spec$var),
    stop(sprintf("Chi-squared does not support specification type '%s'", cls))
  )
}
