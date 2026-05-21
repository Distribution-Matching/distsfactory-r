# Binomial distribution. Parameters: size = n (positive integer), prob = p.
# mean = n*p, var = n*p*(1-p).

binomial_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("Binomial mean and var must be positive")
  if (var >= mean) stop("Binomial requires var < mean")
  p <- 1 - var / mean
  if (!(0 < p && p < 1)) stop("Binomial mean+var implies p out of (0, 1)")
  n_real <- mean / p
  n <- round(n_real)
  if (abs(n - n_real) > 1e-6 || n < 1)
    stop(sprintf("Binomial mean+var requires integer n; got %g", n_real))
  new_dist("binomial", list(size = n, prob = p), dbinom, pbinom, qbinom, rbinom)
}

binomial_exists_mean_var <- function(mean, var) {
  if (!is.finite(mean) || !is.finite(var) || mean <= 0 || var <= 0) return(FALSE)
  if (var >= mean) return(FALSE)
  p <- 1 - var / mean
  if (!(0 < p && p < 1)) return(FALSE)
  n_real <- mean / p
  n_real >= 1 - 1e-9 && abs(n_real - round(n_real)) <= 1e-6
}

binomial_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = binomial_from_mean_var(spec$mean, spec$var),
    stop(sprintf("Binomial does not support specification type '%s'", cls))
  )
}
