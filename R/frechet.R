# Frechet (Type II extreme value) distribution. Parameters: shape alpha > 0,
# scale sigma > 0. Support x >= 0.
# pdf: (alpha/sigma) * (sigma/x)^(alpha+1) * exp(-(sigma/x)^alpha)   for x > 0
# cdf: exp(-(sigma/x)^alpha)
# qf:  sigma * (-log(p))^(-1/alpha)
# mean = sigma * gamma(1 - 1/alpha)         for alpha > 1
# var  = sigma^2 * (gamma(1-2/alpha) - gamma(1-1/alpha)^2)  for alpha > 2

dfrechet_ <- function(x, shape = 1, scale = 1, log = FALSE) {
  if (shape <= 0 || scale <= 0) stop("dfrechet: shape and scale must be positive")
  out <- ifelse(x <= 0, -Inf,
                log(shape) - log(scale) + (shape + 1) * log(scale / x) - (scale / x)^shape)
  if (log) out else exp(out)
}

pfrechet_ <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (shape <= 0 || scale <= 0) stop("pfrechet: shape and scale must be positive")
  p <- ifelse(q <= 0, 0, exp(-(scale / q)^shape))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qfrechet_ <- function(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (shape <= 0 || scale <= 0) stop("qfrechet: shape and scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  scale * (-log(p))^(-1 / shape)
}

rfrechet_ <- function(n, shape = 1, scale = 1) {
  qfrechet_(runif(n), shape = shape, scale = scale)
}

frechet_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("Frechet mean and var must be positive")
  cv2 <- var / mean^2
  # gamma(1-2/alpha)/gamma(1-1/alpha)^2 - 1 = cv2; requires alpha > 2.
  f <- function(log_a) {
    alpha <- 2 + exp(log_a)
    g1 <- gamma(1 - 1 / alpha)
    g2 <- gamma(1 - 2 / alpha)
    g2 / g1^2 - 1 - cv2
  }
  log_a <- find_root_1d(f, x0 = 0)
  alpha <- 2 + exp(log_a)
  sigma <- mean / gamma(1 - 1 / alpha)
  new_dist("frechet", list(shape = alpha, scale = sigma),
           dfrechet_, pfrechet_, qfrechet_, rfrechet_)
}

frechet_from_two_quantiles <- function(p1, q1, p2, q2) {
  if (q1 <= 0 || q2 <= 0) stop("Frechet quantile targets must be positive")
  # log(q_p / sigma) = -(1/alpha) * log(-log(p))
  # log(q1/q2) = (1/alpha) * (log(-log p2) - log(-log p1))   (sign: p<1 => -log(p)>0)
  L1 <- log(-log(p1)); L2 <- log(-log(p2))
  if (isTRUE(all.equal(L1, L2))) stop("Frechet two-quantile: probabilities must differ")
  inv_a <- log(q1 / q2) / (L2 - L1)
  if (inv_a <= 0) stop("Frechet two-quantile: implied shape must be positive")
  alpha <- 1 / inv_a
  sigma <- q1 * (-log(p1))^(1 / alpha)
  new_dist("frechet", list(shape = alpha, scale = sigma),
           dfrechet_, pfrechet_, qfrechet_, rfrechet_)
}

frechet_exists_mean_var <- function(mean, var) {
  if (!is.finite(mean) || !is.finite(var) || mean <= 0 || var <= 0) return(FALSE)
  # var > 0 with mean finite already requires alpha > 2; the cv2 root search
  # determines feasibility numerically. Accept and let the constructor raise
  # if the root finder fails.
  TRUE
}

frechet_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec     = frechet_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec = frechet_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    stop(sprintf("Frechet does not support specification type '%s'", cls))
  )
}
