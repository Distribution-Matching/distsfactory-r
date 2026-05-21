# Negative binomial. Base R: dnbinom(x, size = r, prob = p).
# Distribution counts failures before the r-th success.
# mean = r*(1-p)/p; var = r*(1-p)/p^2.
# So: 1/p = var/mean -> p = mean/var; then r = mean*p/(1-p).

nbinom_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("NegativeBinomial mean and var must be positive")
  if (var <= mean) stop("NegativeBinomial requires var > mean")
  p <- mean / var
  if (!(0 < p && p < 1)) stop("NegativeBinomial mean+var implies p out of (0, 1)")
  r <- mean * p / (1 - p)
  new_dist("negative_binomial", list(size = r, prob = p),
           dnbinom, pnbinom, qnbinom, rnbinom)
}

nbinom_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && mean > 0 && var > mean
}

nbinom_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = nbinom_from_mean_var(spec$mean, spec$var),
    stop(sprintf("NegativeBinomial does not support specification type '%s'", cls))
  )
}
