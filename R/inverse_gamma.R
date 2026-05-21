# Inverse Gamma distribution. If Y ~ Gamma(shape=alpha, scale=theta),
# then X = theta / Y ~ InvGamma(alpha, theta) in the (shape, scale) form, but
# we use the standard (shape, scale) parameterization where the pdf is
#   f(x) = beta^alpha / Gamma(alpha) * x^(-alpha-1) * exp(-beta/x),  x > 0
# with shape = alpha > 0 and scale = beta > 0.
# mean   = beta / (alpha - 1)            for alpha > 1
# var    = beta^2 / ((alpha-1)^2 * (alpha-2))  for alpha > 2
# mode   = beta / (alpha + 1)
# qf(p)  = beta / qgamma(1-p, shape=alpha)  (since X = beta/Y, Y~Gamma(alpha,1))

dinvgamma_ <- function(x, shape = 1, scale = 1, log = FALSE) {
  if (shape <= 0 || scale <= 0) stop("dinvgamma: shape and scale must be positive")
  out <- ifelse(x <= 0, -Inf,
                shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x)
  if (log) out else exp(out)
}

pinvgamma_ <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (shape <= 0 || scale <= 0) stop("pinvgamma: shape and scale must be positive")
  # X = scale/Y, Y ~ Gamma(shape, 1); X > q  iff  Y < scale/q
  p <- ifelse(q <= 0, 0, 1 - pgamma(scale / q, shape = shape, scale = 1))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qinvgamma_ <- function(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (shape <= 0 || scale <= 0) stop("qinvgamma: shape and scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  scale / qgamma(1 - p, shape = shape, scale = 1)
}

rinvgamma_ <- function(n, shape = 1, scale = 1) {
  scale / rgamma(n, shape = shape, scale = 1)
}

invgamma_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("InverseGamma mean and var must be positive")
  # var / mean^2 = 1 / (alpha - 2)
  R <- var / mean^2
  alpha <- 1 / R + 2
  if (alpha <= 2) stop("InverseGamma mean+var infeasible")
  beta <- mean * (alpha - 1)
  new_dist("inverse_gamma", list(shape = alpha, scale = beta),
           dinvgamma_, pinvgamma_, qinvgamma_, rinvgamma_)
}

invgamma_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && mean > 0 && var > 0
}

invgamma_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = invgamma_from_mean_var(spec$mean, spec$var),
    stop(sprintf("InverseGamma does not support specification type '%s'", cls))
  )
}
