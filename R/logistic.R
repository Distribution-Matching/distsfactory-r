# Logistic distribution matching
# Base R parameterization: location, scale
# mean = location, var = scale^2 * pi^2 / 3

logistic_from_mean_var <- function(mu, var) {
  s <- sqrt(3 * var / pi^2)
  new_dist("logistic", list(location = mu, scale = s), dlogis, plogis, qlogis, rlogis)
}

logistic_from_two_quantiles <- function(p1, q1, p2, q2) {
  # location-scale: quantile(p) = location + scale * log(p/(1-p))
  z1 <- log(p1 / (1 - p1))
  z2 <- log(p2 / (1 - p2))
  s <- (q2 - q1) / (z2 - z1)
  mu <- q1 - s * z1
  new_dist("logistic", list(location = mu, scale = s), dlogis, plogis, qlogis, rlogis)
}

logistic_from_mode_iqr <- function(mode, iqr) {
  # mode = location, IQR = 2 * s * ln(3)
  s <- iqr / (2 * log(3))
  new_dist("logistic", list(location = mode, scale = s), dlogis, plogis, qlogis, rlogis)
}

logistic_from_mean_quantile <- function(mu, p, q) {
  z <- log(p / (1 - p))
  if (abs(z) < 1e-12) {
    # p = 0.5: quantile is median, which equals mean for Logistic
    if (abs(mu - q) / max(abs(mu), 1e-12) <= 1e-6) {
      stop("Logistic mean and median are always equal; need an additional constraint to determine scale")
    } else {
      stop(sprintf("Logistic mean must equal median, but got mean=%g, median=%g", mu, q))
    }
  }
  s <- (q - mu) / z
  if (s <= 0) {
    stop(sprintf(
      "Cannot construct Logistic with mean=%g and quantile(%g)=%g: scale would be non-positive",
      mu, p, q
    ))
  }
  new_dist("logistic", list(location = mu, scale = s), dlogis, plogis, qlogis, rlogis)
}

logistic_exists_mean_var <- function(mu, var) {
  var > 0
}

logistic_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = logistic_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec = logistic_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    ModeIQRSpec = logistic_from_mode_iqr(spec$mode, spec$iqr),
    MeanQuantileSpec = logistic_from_mean_quantile(spec$mean, spec$p, spec$q),
    stop(sprintf("Logistic does not support specification type '%s'", cls))
  )
}
