# Beta distribution matching
# Base R parameterization: shape1 (alpha), shape2 (beta) on [0, 1]
# mean = alpha / (alpha + beta)
# var = alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))

beta_from_mean_var <- function(mu, var) {
  S <- mu * (1 - mu) / var - 1
  if (S <= 0) {
    stop(sprintf("Beta variance too large: var=%g >= mu*(1-mu)=%g", var, mu * (1 - mu)))
  }
  alpha <- mu * S
  beta_param <- (1 - mu) * S
  new_dist("beta", list(shape1 = alpha, shape2 = beta_param), dbeta, pbeta, qbeta, rbeta)
}

beta_from_mean_mode <- function(mu, mode) {
  if (abs(mu - mode) < 1e-12) {
    stop("Beta mean and mode cannot be equal (symmetric case is underdetermined)")
  }
  alpha <- mu * (2 * mode - 1) / (mode - mu)
  beta_param <- alpha * (1 - mu) / mu
  if (alpha <= 1 || beta_param <= 1) {
    stop(sprintf(
      "Beta from mean=%g, mode=%g gives alpha=%.4f, beta=%.4f; both must be > 1 for mode to exist",
      mu, mode, alpha, beta_param
    ))
  }
  new_dist("beta", list(shape1 = alpha, shape2 = beta_param), dbeta, pbeta, qbeta, rbeta)
}

beta_from_two_quantiles <- function(p1, q1, p2, q2) {
  # Initial guess via normal approximation
  z1 <- qnorm(p1)
  z2 <- qnorm(p2)
  mu_est <- (q1 * z2 - q2 * z1) / (z2 - z1)
  sigma_est <- (q2 - q1) / (z2 - z1)
  var_est <- sigma_est^2

  # Clamp to valid Beta range
  mu_est <- max(0.01, min(0.99, mu_est))
  var_est <- min(var_est, mu_est * (1 - mu_est) * 0.9)
  S0 <- mu_est * (1 - mu_est) / var_est - 1
  alpha0 <- max(0.5, mu_est * S0)
  beta0 <- max(0.5, (1 - mu_est) * S0)

  F_resid <- function(x) {
    a <- exp(x[1])
    b <- exp(x[2])
    c(qbeta(p1, a, b) - q1, qbeta(p2, a, b) - q2)
  }

  x <- newton_2d(F_resid, c(log(alpha0), log(beta0)))
  alpha <- exp(x[1])
  beta_param <- exp(x[2])
  new_dist("beta", list(shape1 = alpha, shape2 = beta_param), dbeta, pbeta, qbeta, rbeta)
}

beta_from_mean_quantile <- function(mu, p, q) {
  # S = alpha + beta; alpha = mu*S, beta = (1-mu)*S
  f <- function(log_S) {
    S <- exp(log_S)
    a <- mu * S
    b <- (1 - mu) * S
    qbeta(p, a, b) - q
  }
  log_S <- find_root_1d(f, x0 = 1)
  S <- exp(log_S)
  alpha <- mu * S
  beta_param <- (1 - mu) * S
  new_dist("beta", list(shape1 = alpha, shape2 = beta_param), dbeta, pbeta, qbeta, rbeta)
}

beta_exists_mean_var <- function(mu, var) {
  mu > 0 && mu < 1 && var > 0 && var < mu * (1 - mu)
}

beta_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = beta_from_mean_var(spec$mean, spec$var),
    MeanModeSpec = beta_from_mean_mode(spec$mean, spec$mode),
    TwoQuantileSpec = beta_from_two_quantiles(spec$p1, spec$q1, spec$p2, spec$q2),
    MeanQuantileSpec = beta_from_mean_quantile(spec$mean, spec$p, spec$q),
    stop(sprintf("Beta does not support specification type '%s'", cls))
  )
}
