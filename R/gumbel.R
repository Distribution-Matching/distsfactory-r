# Gumbel (right-skewed extreme value) distribution.
# Parameters: location = mu, scale = beta > 0.
# pdf: (1/beta) * exp(-(z + exp(-z)))   where z = (x - mu)/beta
# cdf: exp(-exp(-z))
# qf:  mu - beta * log(-log(p))
# mean   = mu + beta * EulerGamma
# var    = pi^2 * beta^2 / 6
# median = mu - beta * log(log(2))
# mode   = mu

.EULER_MASCHERONI <- -digamma(1)  # ~ 0.5772156649

dgumbel_ <- function(x, location = 0, scale = 1, log = FALSE) {
  if (scale <= 0) stop("dgumbel: scale must be positive")
  z <- (x - location) / scale
  out <- -log(scale) - z - exp(-z)
  if (log) out else exp(out)
}

pgumbel_ <- function(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("pgumbel: scale must be positive")
  z <- (q - location) / scale
  p <- exp(-exp(-z))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qgumbel_ <- function(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("qgumbel: scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  location - scale * log(-log(p))
}

rgumbel_ <- function(n, location = 0, scale = 1) {
  u <- runif(n)
  location - scale * log(-log(u))
}

# Quantile linear coefficient: q_p = mu + beta * c(p) where c(p) = -log(-log(p))
.gumbel_c <- function(p) -log(-log(p))

gumbel_from_mean_var <- function(mean, var) {
  if (var <= 0) stop("Gumbel variance must be positive")
  beta <- sqrt(6 * var) / pi
  mu <- mean - beta * .EULER_MASCHERONI
  new_dist("gumbel", list(location = mu, scale = beta),
           dgumbel_, pgumbel_, qgumbel_, rgumbel_)
}

gumbel_from_two_quantiles <- function(p1, q1, p2, q2) {
  c1 <- .gumbel_c(p1); c2 <- .gumbel_c(p2)
  if (isTRUE(all.equal(c1, c2)))
    stop("Gumbel two-quantile: probabilities must differ")
  beta <- (q2 - q1) / (c2 - c1)
  if (beta <= 0) stop("Gumbel two-quantile: implied scale must be positive")
  mu <- q1 - beta * c1
  new_dist("gumbel", list(location = mu, scale = beta),
           dgumbel_, pgumbel_, qgumbel_, rgumbel_)
}

gumbel_from_mean_quantile <- function(mean, p, q) {
  c_p <- .gumbel_c(p)
  # q = mu + beta*c_p; mean = mu + beta*EulerGamma
  # q - mean = beta * (c_p - EulerGamma)
  denom <- c_p - .EULER_MASCHERONI
  if (isTRUE(all.equal(denom, 0)))
    stop("Gumbel mean+quantile: probability coincides with mean")
  beta <- (q - mean) / denom
  if (beta <= 0) stop("Gumbel mean+quantile: implied scale must be positive")
  mu <- mean - beta * .EULER_MASCHERONI
  new_dist("gumbel", list(location = mu, scale = beta),
           dgumbel_, pgumbel_, qgumbel_, rgumbel_)
}

gumbel_from_mean_mode <- function(mean, mode) {
  beta <- (mean - mode) / .EULER_MASCHERONI
  if (beta <= 0) stop("Gumbel: mean must exceed mode")
  new_dist("gumbel", list(location = mode, scale = beta),
           dgumbel_, pgumbel_, qgumbel_, rgumbel_)
}

gumbel_from_mode_var <- function(mode, var) {
  if (var <= 0) stop("Gumbel variance must be positive")
  beta <- sqrt(6 * var) / pi
  new_dist("gumbel", list(location = mode, scale = beta),
           dgumbel_, pgumbel_, qgumbel_, rgumbel_)
}

gumbel_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && var > 0
}

gumbel_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec      = gumbel_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec  = gumbel_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    MeanQuantileSpec = gumbel_from_mean_quantile(spec$mean, spec$p, spec$q),
    MeanModeSpec     = gumbel_from_mean_mode(spec$mean, spec$mode),
    ModeVarSpec      = gumbel_from_mode_var(spec$mode, spec$var),
    stop(sprintf("Gumbel does not support specification type '%s'", cls))
  )
}
