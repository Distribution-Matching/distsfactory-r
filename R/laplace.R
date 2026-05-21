# Laplace distribution (double exponential), location-scale.
# Parameters: location = mu, scale = b > 0.
# pdf:   (1/(2b)) * exp(-|x-mu|/b)
# cdf:   0.5 + 0.5 * sign(x-mu) * (1 - exp(-|x-mu|/b))
# qf:    mu + b * psi(p) where psi(p) = sign(p-0.5)*(-log(1 - 2*|p-0.5|))
#         = log(2p) for p < 0.5; = -log(2(1-p)) for p > 0.5; 0 at 0.5
# mean = mu, var = 2*b^2, median = mode = mu.

dlaplace_ <- function(x, location = 0, scale = 1, log = FALSE) {
  if (scale <= 0) stop("dlaplace: scale must be positive")
  z <- abs(x - location) / scale
  out <- -log(2 * scale) - z
  if (log) out else exp(out)
}

plaplace_ <- function(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("plaplace: scale must be positive")
  z <- (q - location) / scale
  p <- ifelse(z < 0,
              0.5 * exp(z),
              1 - 0.5 * exp(-z))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qlaplace_ <- function(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("qlaplace: scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  ifelse(p < 0.5,
         location + scale * log(2 * p),
         location - scale * log(2 * (1 - p)))
}

rlaplace_ <- function(n, location = 0, scale = 1) {
  u <- runif(n) - 0.5
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# Quantile-function linear coefficient: q_p = mu + b * psi(p).
.laplace_psi <- function(p) {
  ifelse(p < 0.5, log(2 * p), -log(2 * (1 - p)))
}

laplace_from_mean_var <- function(mean, var) {
  if (var <= 0) stop("Laplace variance must be positive")
  new_dist("laplace", list(location = mean, scale = sqrt(var / 2)),
           dlaplace_, plaplace_, qlaplace_, rlaplace_)
}

laplace_from_two_quantiles <- function(p1, q1, p2, q2) {
  psi1 <- .laplace_psi(p1); psi2 <- .laplace_psi(p2)
  if (isTRUE(all.equal(psi1, psi2)))
    stop("Laplace two-quantile: probabilities must differ enough to pin scale")
  b <- (q2 - q1) / (psi2 - psi1)
  if (b <= 0) stop("Laplace two-quantile: implied scale must be positive")
  mu <- q1 - b * psi1
  new_dist("laplace", list(location = mu, scale = b),
           dlaplace_, plaplace_, qlaplace_, rlaplace_)
}

laplace_from_mean_quantile <- function(mean, p, q) {
  psi <- .laplace_psi(p)
  if (isTRUE(all.equal(psi, 0)))
    stop("Laplace mean+median underdetermines scale")
  b <- (q - mean) / psi
  if (b <= 0) stop("Laplace mean+quantile: implied scale must be positive")
  new_dist("laplace", list(location = mean, scale = b),
           dlaplace_, plaplace_, qlaplace_, rlaplace_)
}

laplace_from_mode_var <- function(mode, var) laplace_from_mean_var(mode, var)

laplace_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && var > 0
}

laplace_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec      = laplace_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec  = laplace_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    MeanQuantileSpec = laplace_from_mean_quantile(spec$mean, spec$p, spec$q),
    ModeVarSpec      = laplace_from_mode_var(spec$mode, spec$var),
    stop(sprintf("Laplace does not support specification type '%s'", cls))
  )
}
