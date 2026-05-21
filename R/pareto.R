# Pareto (Type I) distribution. Parameters: shape alpha > 0, scale xm > 0.
# Support: x >= xm.
# pdf: alpha * xm^alpha / x^(alpha+1)
# cdf: 1 - (xm/x)^alpha
# qf:  xm * (1-p)^(-1/alpha)
# mean = alpha*xm / (alpha-1)        for alpha > 1
# var  = xm^2 * alpha / ((alpha-1)^2 * (alpha-2))  for alpha > 2

dpareto_ <- function(x, shape = 1, scale = 1, log = FALSE) {
  if (shape <= 0 || scale <= 0) stop("dpareto: shape and scale must be positive")
  out <- ifelse(x < scale, -Inf,
                log(shape) + shape * log(scale) - (shape + 1) * log(x))
  if (log) out else exp(out)
}

ppareto_ <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (shape <= 0 || scale <= 0) stop("ppareto: shape and scale must be positive")
  p <- ifelse(q < scale, 0, 1 - (scale / q)^shape)
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qpareto_ <- function(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (shape <= 0 || scale <= 0) stop("qpareto: shape and scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  scale * (1 - p)^(-1 / shape)
}

rpareto_ <- function(n, shape = 1, scale = 1) {
  qpareto_(runif(n), shape = shape, scale = scale)
}

pareto_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("Pareto mean and var must be positive")
  # mu = alpha*xm/(alpha-1); var = xm^2 * alpha / ((alpha-1)^2 * (alpha-2))
  # Ratio var/mu^2 = (alpha-1)^2 * (alpha-2) * alpha / (alpha^2 * (alpha-1)^2 * (alpha-2))
  # Actually:  mu^2 = alpha^2 * xm^2 / (alpha-1)^2
  # var / mu^2 = xm^2*alpha / ((alpha-1)^2*(alpha-2)) * (alpha-1)^2 / (alpha^2 * xm^2)
  #            = 1 / (alpha * (alpha-2))
  # So alpha*(alpha-2) = 1/(var/mu^2)  ->  alpha^2 - 2*alpha - mu^2/var = 0
  R <- mean^2 / var
  disc <- 4 + 4 * R
  alpha <- 1 + sqrt(1 + R)  # take positive root, requires alpha > 2
  if (alpha <= 2) stop("Pareto mean+var: implied shape must exceed 2")
  xm <- mean * (alpha - 1) / alpha
  new_dist("pareto", list(shape = alpha, scale = xm),
           dpareto_, ppareto_, qpareto_, rpareto_)
}

pareto_from_two_quantiles <- function(p1, q1, p2, q2) {
  if (q1 <= 0 || q2 <= 0) stop("Pareto quantile targets must be positive")
  # q_p = xm * (1-p)^(-1/alpha)
  # log(q1/q2) = (1/alpha) * log((1-p2)/(1-p1))
  num <- log(q1 / q2)
  den <- log((1 - p2) / (1 - p1))
  if (isTRUE(all.equal(den, 0))) stop("Pareto two-quantile: probabilities must differ")
  inv_alpha <- num / den
  if (inv_alpha <= 0) stop("Pareto two-quantile: implied shape must be positive")
  alpha <- 1 / inv_alpha
  xm <- q1 * (1 - p1)^(1 / alpha)
  new_dist("pareto", list(shape = alpha, scale = xm),
           dpareto_, ppareto_, qpareto_, rpareto_)
}

pareto_exists_mean_var <- function(mean, var) {
  if (!is.finite(mean) || !is.finite(var) || mean <= 0 || var <= 0) return(FALSE)
  R <- mean^2 / var
  alpha <- 1 + sqrt(1 + R)
  alpha > 2
}

pareto_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec     = pareto_from_mean_var(spec$mean, spec$var),
    TwoQuantileSpec = pareto_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    stop(sprintf("Pareto does not support specification type '%s'", cls))
  )
}
