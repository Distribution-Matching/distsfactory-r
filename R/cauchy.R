# Cauchy distribution. Two parameters: location (median = mode) and scale > 0.
# Base R parameterization: dcauchy(x, location, scale).
# pdf:    1 / (pi * scale * (1 + ((x - location)/scale)^2))
# cdf:    1/2 + atan((x - location)/scale) / pi
# qf(p):  location + scale * tan(pi * (p - 1/2))
# median = mode = location;  IQR = 2 * scale.
# Mean and variance are *undefined* — construction is via quantiles only.

cauchy_from_two_quantiles <- function(p1, q1, p2, q2) {
  # qcauchy(p) = location + scale * tan(pi*(p - 1/2)). Linear in (loc, scale).
  t1 <- tan(pi * (p1 - 0.5))
  t2 <- tan(pi * (p2 - 0.5))
  if (isTRUE(all.equal(t1, t2)))
    stop("Cauchy two-quantile: probabilities must differ")
  scale <- (q2 - q1) / (t2 - t1)
  if (scale <= 0) stop("Cauchy: implied scale must be positive")
  location <- q1 - scale * t1
  new_dist("cauchy", list(location = location, scale = scale),
           dcauchy, pcauchy, qcauchy, rcauchy)
}

cauchy_from_median_iqr <- function(median, iqr) {
  if (iqr <= 0) stop("Cauchy IQR must be positive")
  # IQR = q(0.75) - q(0.25) = scale * (tan(pi/4) - tan(-pi/4)) = 2 * scale
  new_dist("cauchy", list(location = median, scale = iqr / 2),
           dcauchy, pcauchy, qcauchy, rcauchy)
}

# Cauchy has no finite moments — `dist_exists("cauchy", mean=, var=)` is
# always FALSE, matching Julia's `Cauchy: no finite mean or variance`.
cauchy_exists_mean_var <- function(mean, var) FALSE

cauchy_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    TwoQuantileSpec = cauchy_from_two_quantiles(spec$p1, spec$q1,
                                                spec$p2, spec$q2),
    MeanVarSpec = stop(
      "Cauchy has no finite mean or variance; construct via quantiles instead"
    ),
    MeanSpec = stop(
      "Cauchy has no finite mean; construct via quantiles instead"
    ),
    VarSpec = stop(
      "Cauchy has no finite variance; construct via quantiles instead"
    ),
    stop(sprintf("Cauchy does not support specification type '%s'", cls))
  )
}
