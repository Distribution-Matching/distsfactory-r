# Chi distribution. X = sqrt(Y) where Y ~ Chisq(k); one parameter df = k > 0.
# pdf(x) = x^(k-1) * exp(-x^2/2) / (2^(k/2-1) * Gamma(k/2))   for x > 0
# E[X]  = sqrt(2) * Gamma((k+1)/2) / Gamma(k/2)
# Var[X] = k - E[X]^2
# cdf(x) = P(Y <= x^2) = pchisq(x^2, df = k)
# qf(p) = sqrt(qchisq(p, df = k))

dchi_ <- function(x, df = 1, log = FALSE) {
  if (df <= 0) stop("dchi: df must be positive")
  out <- ifelse(x <= 0, -Inf,
                (df - 1) * log(x) - x^2 / 2 - (df / 2 - 1) * log(2) - lgamma(df / 2))
  if (log) out else exp(out)
}

pchi_ <- function(q, df = 1, lower.tail = TRUE, log.p = FALSE) {
  if (df <= 0) stop("pchi: df must be positive")
  p <- ifelse(q <= 0, 0, pchisq(q^2, df = df))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qchi_ <- function(p, df = 1, lower.tail = TRUE, log.p = FALSE) {
  if (df <= 0) stop("qchi: df must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  sqrt(qchisq(p, df = df))
}

rchi_ <- function(n, df = 1) sqrt(rchisq(n, df = df))

chi_mean_of_df <- function(df) sqrt(2) * exp(lgamma((df + 1) / 2) - lgamma(df / 2))

chi_from_mean <- function(mean) {
  if (mean <= 0) stop("Chi mean must be positive")
  f <- function(log_df) chi_mean_of_df(exp(log_df)) - mean
  log_df <- find_root_1d(f, x0 = 0)
  df <- exp(log_df)
  new_dist("chi", list(df = df), dchi_, pchi_, qchi_, rchi_)
}

chi_from_mean_var <- function(mean, var) {
  if (mean <= 0 || var <= 0) stop("Chi mean and var must be positive")
  # var = df - mean^2, so df = mean^2 + var. Then verify against chi_mean_of_df.
  df <- mean^2 + var
  expected <- chi_mean_of_df(df)
  if (!isTRUE(all.equal(expected, mean, tolerance = 1e-5)))
    stop(sprintf("Chi mean+var infeasible: implied df=%g gives mean=%g, requested %g",
                 df, expected, mean))
  new_dist("chi", list(df = df), dchi_, pchi_, qchi_, rchi_)
}

chi_exists_mean_var <- function(mean, var) {
  if (!is.finite(mean) || !is.finite(var) || mean <= 0 || var <= 0) return(FALSE)
  df <- mean^2 + var
  isTRUE(all.equal(chi_mean_of_df(df), mean, tolerance = 1e-5))
}

chi_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanSpec    = chi_from_mean(spec$mean),
    MeanVarSpec = chi_from_mean_var(spec$mean, spec$var),
    stop(sprintf("Chi does not support specification type '%s'", cls))
  )
}
