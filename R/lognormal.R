# Lognormal distribution matching.
# Base R parameterization: dlnorm(x, meanlog = mu, sdlog = sigma)
# X ~ Lognormal(mu, sigma)  iff  log(X) ~ Normal(mu, sigma)
# mean = exp(mu + sigma^2/2)
# var  = (exp(sigma^2) - 1) * exp(2*mu + sigma^2)
# median = exp(mu)
# mode   = exp(mu - sigma^2)

lognormal_from_mean_var <- function(mean, var) {
  if (mean <= 0) stop(sprintf("Lognormal mean must be positive, got %g", mean))
  if (var <= 0) stop(sprintf("Lognormal variance must be positive, got %g", var))
  # 1 + var/mean^2 = exp(sigma^2)  ->  sigma^2 = log(1 + var/mean^2)
  sigma2 <- log1p(var / mean^2)
  sigma <- sqrt(sigma2)
  mu <- log(mean) - sigma2 / 2
  new_dist("lognormal", list(meanlog = mu, sdlog = sigma),
           dlnorm, plnorm, qlnorm, rlnorm)
}

lognormal_from_two_quantiles <- function(p1, q1, p2, q2) {
  if (q1 <= 0 || q2 <= 0) stop("Lognormal quantile targets must be positive")
  # log of quantiles are normal quantiles
  z1 <- qnorm(p1); z2 <- qnorm(p2)
  if (z2 == z1) stop("two-quantile lognormal: distinct probabilities required")
  sigma <- (log(q2) - log(q1)) / (z2 - z1)
  if (sigma <= 0) stop("lognormal: implied sigma must be positive")
  mu <- log(q1) - sigma * z1
  new_dist("lognormal", list(meanlog = mu, sdlog = sigma),
           dlnorm, plnorm, qlnorm, rlnorm)
}

lognormal_from_mean_quantile <- function(mean, p, q) {
  if (mean <= 0) stop("Lognormal mean must be positive")
  if (q <= 0) stop("Lognormal quantile target must be positive")
  # Solve for sigma > 0 given mean = exp(mu + sigma^2/2) and qnorm(p)*sigma + mu = log(q)
  # From the first: mu = log(mean) - sigma^2/2
  # Substituted:    log(mean) - sigma^2/2 + qnorm(p)*sigma = log(q)
  z <- qnorm(p)
  rhs <- log(q) - log(mean)  # = -sigma^2/2 + z*sigma
  f <- function(sigma) -sigma^2 / 2 + z * sigma - rhs
  # f(0) = -rhs. For p<0.5 and q<mean we expect a root sigma > 0.
  # df/dsigma = -sigma + z = 0 at sigma = z (only positive if z>0).
  # Use a wide-bracket search.
  root <- find_root_1d(f, x0 = 1)
  sigma <- root
  if (sigma <= 0) stop("lognormal: no positive sigma satisfies the given mean+quantile")
  mu <- log(mean) - sigma^2 / 2
  new_dist("lognormal", list(meanlog = mu, sdlog = sigma),
           dlnorm, plnorm, qlnorm, rlnorm)
}

lognormal_from_mean_mode <- function(mean, mode) {
  if (mean <= mode) stop(sprintf("Lognormal mean (%g) must exceed mode (%g)", mean, mode))
  # mean/mode = exp(3*sigma^2/2)  ->  sigma^2 = (2/3) * log(mean/mode)
  sigma2 <- (2 / 3) * log(mean / mode)
  if (sigma2 <= 0) stop("Lognormal mean+mode: implied sigma must be positive")
  sigma <- sqrt(sigma2)
  mu <- log(mean) - sigma2 / 2
  new_dist("lognormal", list(meanlog = mu, sdlog = sigma),
           dlnorm, plnorm, qlnorm, rlnorm)
}

lognormal_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && mean > 0 && var > 0
}

lognormal_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec      = lognormal_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec  = lognormal_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    MeanQuantileSpec = lognormal_from_mean_quantile(spec$mean, spec$p, spec$q),
    MeanModeSpec     = lognormal_from_mean_mode(spec$mean, spec$mode),
    stop(sprintf("Lognormal does not support specification type '%s'", cls))
  )
}
