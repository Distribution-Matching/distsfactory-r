# Rayleigh distribution. One-parameter: scale sigma > 0.
# pdf: (x/sigma^2) * exp(-x^2/(2*sigma^2)),  x >= 0
# cdf: 1 - exp(-x^2/(2*sigma^2))
# qf:  sigma * sqrt(-2*log(1-p))
# mean   = sigma * sqrt(pi/2)
# var    = sigma^2 * (4 - pi)/2
# median = sigma * sqrt(2*log(2))
# mode   = sigma

drayleigh_ <- function(x, scale = 1, log = FALSE) {
  if (scale <= 0) stop("drayleigh: scale must be positive")
  out <- ifelse(x < 0, -Inf,
                log(x) - 2 * log(scale) - x^2 / (2 * scale^2))
  if (log) out else exp(out)
}

prayleigh_ <- function(q, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("prayleigh: scale must be positive")
  p <- ifelse(q < 0, 0, 1 - exp(-q^2 / (2 * scale^2)))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qrayleigh_ <- function(p, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("qrayleigh: scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  scale * sqrt(-2 * log1p(-p))
}

rrayleigh_ <- function(n, scale = 1) {
  qrayleigh_(runif(n), scale = scale)
}

rayleigh_from_mean <- function(mean) {
  if (mean <= 0) stop("Rayleigh mean must be positive")
  sigma <- mean / sqrt(pi / 2)
  new_dist("rayleigh", list(scale = sigma),
           drayleigh_, prayleigh_, qrayleigh_, rrayleigh_)
}

rayleigh_from_mean_var <- function(mean, var) {
  if (mean <= 0) stop("Rayleigh mean must be positive")
  expected_var <- mean^2 * (4 - pi) / pi
  if (!isTRUE(all.equal(var, expected_var, tolerance = 1e-6)))
    stop(sprintf("Rayleigh requires var = mean^2*(4-pi)/pi = %g; got %g",
                 expected_var, var))
  rayleigh_from_mean(mean)
}

rayleigh_from_mode <- function(mode) {
  if (mode <= 0) stop("Rayleigh mode must be positive")
  new_dist("rayleigh", list(scale = mode),
           drayleigh_, prayleigh_, qrayleigh_, rrayleigh_)
}

rayleigh_from_median <- function(median) {
  if (median <= 0) stop("Rayleigh median must be positive")
  sigma <- median / sqrt(2 * log(2))
  new_dist("rayleigh", list(scale = sigma),
           drayleigh_, prayleigh_, qrayleigh_, rrayleigh_)
}

rayleigh_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && mean > 0 && var > 0 &&
    isTRUE(all.equal(var, mean^2 * (4 - pi) / pi, tolerance = 1e-6))
}

rayleigh_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanSpec     = rayleigh_from_mean(spec$mean),
    MeanVarSpec  = rayleigh_from_mean_var(spec$mean, spec$var),
    ModeSpec     = rayleigh_from_mode(spec$mode),
    QuantileSpec = if (isTRUE(all.equal(spec$p, 0.5))) rayleigh_from_median(spec$q)
                   else stop("Rayleigh single-quantile only supports median"),
    stop(sprintf("Rayleigh does not support specification type '%s'", cls))
  )
}
