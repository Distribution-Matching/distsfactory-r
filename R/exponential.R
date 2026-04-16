# Exponential distribution matching
# Base R parameterization: rate (mean = 1/rate)

exponential_from_mean <- function(mu) {
  if (mu <= 0) stop(sprintf("Exponential mean must be positive, got %g", mu))
  new_dist("exponential", list(rate = 1 / mu), dexp, pexp, qexp, rexp)
}

exponential_from_mean_var <- function(mu, var) {
  expected_var <- mu^2
  if (abs(var - expected_var) / max(expected_var, 1e-12) > 1e-6) {
    stop(sprintf(
      "Exponential variance must equal mean^2. Got mean=%g, var=%g, expected var=%g",
      mu, var, expected_var
    ))
  }
  exponential_from_mean(mu)
}

exponential_from_var <- function(var) {
  mu <- sqrt(var)
  exponential_from_mean(mu)
}

exponential_from_quantile <- function(p, q) {
  theta <- -q / log(1 - p)
  new_dist("exponential", list(rate = 1 / theta), dexp, pexp, qexp, rexp)
}

exponential_exists_mean_var <- function(mu, var) {
  mu > 0 && var > 0 && abs(var - mu^2) / max(mu^2, 1e-12) <= 1e-6
}

exponential_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = exponential_from_mean_var(spec$mean, spec$var),
    MeanSpec = exponential_from_mean(spec$mean),
    VarSpec = exponential_from_var(spec$var),
    QuantileSpec = exponential_from_quantile(spec$p, spec$q),
    stop(sprintf("Exponential does not support specification type '%s'", cls))
  )
}
