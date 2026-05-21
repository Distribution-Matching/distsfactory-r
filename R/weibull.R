# Weibull distribution matching.
# Base R parameterization: dweibull(x, shape, scale)
# Let k = shape, lambda = scale.
# mean = lambda * gamma(1 + 1/k)
# var  = lambda^2 * (gamma(1 + 2/k) - gamma(1 + 1/k)^2)
# qweibull(p) = lambda * (-log(1 - p))^(1/k)
# mode (k > 1) = lambda * ((k-1)/k)^(1/k); mode = 0 for k <= 1.

weibull_from_mean_var <- function(mean, var) {
  if (mean <= 0) stop(sprintf("Weibull mean must be positive, got %g", mean))
  if (var <= 0) stop(sprintf("Weibull variance must be positive, got %g", var))
  cv2 <- var / mean^2
  # Solve gamma(1+2/k) / gamma(1+1/k)^2 - 1 = cv2 for k > 0.
  f <- function(log_k) {
    k <- exp(log_k)
    g1 <- gamma(1 + 1 / k)
    g2 <- gamma(1 + 2 / k)
    g2 / g1^2 - 1 - cv2
  }
  log_k <- find_root_1d(f, x0 = 0)
  k <- exp(log_k)
  lambda <- mean / gamma(1 + 1 / k)
  new_dist("weibull", list(shape = k, scale = lambda),
           dweibull, pweibull, qweibull, rweibull)
}

weibull_from_two_quantiles <- function(p1, q1, p2, q2) {
  if (q1 <= 0 || q2 <= 0) stop("Weibull quantile targets must be positive")
  L1 <- -log1p(-p1)
  L2 <- -log1p(-p2)
  if (L1 <= 0 || L2 <= 0) stop("Weibull two-quantile: invalid probabilities")
  if (isTRUE(all.equal(L1, L2))) stop("Weibull two-quantile: probabilities must differ")
  # 1/k = log(q1/q2) / log(L1/L2)
  inv_k <- log(q1 / q2) / log(L1 / L2)
  if (inv_k <= 0) stop("Weibull two-quantile: implied shape must be positive")
  k <- 1 / inv_k
  lambda <- q1 / L1^inv_k
  new_dist("weibull", list(shape = k, scale = lambda),
           dweibull, pweibull, qweibull, rweibull)
}

weibull_from_mean_quantile <- function(mean, p, q) {
  if (mean <= 0) stop("Weibull mean must be positive")
  if (q <= 0) stop("Weibull quantile target must be positive")
  L <- -log1p(-p)
  if (L <= 0) stop("Weibull mean+quantile: invalid probability")
  ratio <- q / mean
  f <- function(log_k) {
    k <- exp(log_k)
    L^(1 / k) / gamma(1 + 1 / k) - ratio
  }
  log_k <- find_root_1d(f, x0 = 0)
  k <- exp(log_k)
  lambda <- mean / gamma(1 + 1 / k)
  new_dist("weibull", list(shape = k, scale = lambda),
           dweibull, pweibull, qweibull, rweibull)
}

weibull_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && mean > 0 && var > 0
}

weibull_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec      = weibull_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec  = weibull_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    MeanQuantileSpec = weibull_from_mean_quantile(spec$mean, spec$p, spec$q),
    stop(sprintf("Weibull does not support specification type '%s'", cls))
  )
}
