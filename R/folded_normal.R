# Folded Normal. X = |Y| where Y ~ N(mu_underlying, sigma_underlying^2).
# Two parameters: location = mu, scale = sigma > 0. Both refer to the
# underlying Normal *before* folding.
# pdf(x) = (1/sigma) * (phi((x-mu)/sigma) + phi((-x-mu)/sigma)),  x >= 0
# cdf(x) = Phi((x-mu)/sigma) + Phi((x+mu)/sigma) - 1                (since
#         pre-folded N is symmetric about mu, this is the absolute-value cdf)
# mean = sigma*sqrt(2/pi)*exp(-mu^2/(2*sigma^2)) + mu*(1 - 2*Phi(-mu/sigma))
# var  = mu^2 + sigma^2 - mean^2

dfoldednorm_ <- function(x, location = 0, scale = 1, log = FALSE) {
  if (scale <= 0) stop("dfoldednorm: scale must be positive")
  out <- ifelse(x < 0, -Inf,
                log(dnorm((x - location) / scale) + dnorm((-x - location) / scale)) -
                  log(scale))
  if (log) out else exp(out)
}

pfoldednorm_ <- function(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("pfoldednorm: scale must be positive")
  p <- ifelse(q < 0, 0,
              pnorm((q - location) / scale) + pnorm((q + location) / scale) - 1)
  p <- pmin(pmax(p, 0), 1)
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qfoldednorm_ <- function(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("qfoldednorm: scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  out <- numeric(length(p))
  for (i in seq_along(p)) {
    target <- p[i]
    if (target <= 0) { out[i] <- 0; next }
    if (target >= 1) { out[i] <- Inf; next }
    f <- function(q) pfoldednorm_(q, location, scale) - target
    out[i] <- uniroot(f, lower = 0, upper = max(1, abs(location) + 10 * scale),
                      extendInt = "upX", tol = 1e-12)$root
  }
  out
}

rfoldednorm_ <- function(n, location = 0, scale = 1) {
  abs(rnorm(n, mean = location, sd = scale))
}

.foldednorm_mean_var <- function(mu, sigma) {
  m <- sigma * sqrt(2 / pi) * exp(-mu^2 / (2 * sigma^2)) +
       mu * (1 - 2 * pnorm(-mu / sigma))
  v <- mu^2 + sigma^2 - m^2
  c(mean = m, var = v)
}

foldednorm_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("FoldedNormal mean and var must be positive")
  # 2-D root-find on (mu, log_sigma) for the system (mean, var) match.
  obj <- function(params) {
    mu <- params[1]
    sigma <- exp(params[2])
    mv <- .foldednorm_mean_var(mu, sigma)
    c(mv["mean"] - mean, mv["var"] - var)
  }
  # Initial guess: pure half-normal (mu=0) -> mean = sigma*sqrt(2/pi); var = sigma^2*(1-2/pi).
  sigma0 <- sqrt(max(var / (1 - 2 / pi), 1e-6))
  start <- c(0, log(sigma0))
  sol <- newton_2d(obj, start)
  mu <- sol[1]; sigma <- exp(sol[2])
  new_dist("folded_normal", list(location = mu, scale = sigma),
           dfoldednorm_, pfoldednorm_, qfoldednorm_, rfoldednorm_)
}

foldednorm_exists_mean_var <- function(mean, var) {
  if (!is.finite(mean) || !is.finite(var) || mean <= 0 || var <= 0) return(FALSE)
  # Pure half-normal (mu=0) has the largest var/mean^2 = (pi/2 - 1).
  # Any folded normal must have var/mean^2 <= pi/2 - 1.
  upper_bound <- pi / 2 - 1
  var / mean^2 <= upper_bound + 1e-9
}

foldednorm_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = foldednorm_from_mean_var(spec$mean, spec$var),
    stop(sprintf("FoldedNormal does not support specification type '%s'", cls))
  )
}
