# Closed-form mean and variance for each family, parameterized by the
# distsfactory_dist `$params` list. Used by the partial_dist solver.

dist_mean_from_params <- function(name, params) {
  switch(name,
    gamma       = params$shape / params$rate,
    exponential = 1 / params$rate,
    logistic    = params$location,
    beta        = params$shape1 / (params$shape1 + params$shape2),
    normal      = params$mean,
    lognormal   = exp(params$meanlog + params$sdlog^2 / 2),
    uniform     = (params$min + params$max) / 2,
    weibull     = params$scale * gamma(1 + 1 / params$shape),
    tdist       = 0,
    chisq       = params$df,
    fdist       = params$df2 / (params$df2 - 2),
    laplace     = params$location,
    gumbel      = params$location + params$scale * (-digamma(1)),
    rayleigh    = params$scale * sqrt(pi / 2),
    pareto      = params$shape * params$scale / (params$shape - 1),
    frechet     = params$scale * gamma(1 - 1 / params$shape),
    inverse_gamma = params$scale / (params$shape - 1),
    chi         = sqrt(2) * exp(lgamma((params$df + 1) / 2) - lgamma(params$df / 2)),
    folded_normal = {
      mu <- params$location; sigma <- params$scale
      sigma * sqrt(2 / pi) * exp(-mu^2 / (2 * sigma^2)) +
        mu * (1 - 2 * pnorm(-mu / sigma))
    },
    erlang      = params$shape * params$scale,
    sym_triangular = params$location,
    triangular  = (params$a + params$b + params$c) / 3,
    binomial    = params$size * params$prob,
    poisson     = params$lambda,
    negative_binomial = params$size * (1 - params$prob) / params$prob,
    geometric   = (1 - params$prob) / params$prob,
    discrete_uniform = (params$min + params$max) / 2,
    discrete_sym_triangular = params$mu,
    discrete_triangular = {
      a <- params$a; b <- params$b; c <- params$c
      Z <- (b - a + 2) / 2
      ks <- a:b
      pmf <- ifelse(ks <= c, (ks - a + 1) / (c - a + 1) / Z,
                              (b - ks + 1) / (b - c + 1) / Z)
      sum(ks * pmf)
    },
    stop(sprintf("no closed-form mean wired for %s", name))
  )
}

dist_var_from_params <- function(name, params) {
  switch(name,
    gamma       = params$shape / params$rate^2,
    exponential = 1 / params$rate^2,
    logistic    = (pi^2 / 3) * params$scale^2,
    beta        = {
      a <- params$shape1; b <- params$shape2
      a * b / ((a + b)^2 * (a + b + 1))
    },
    normal      = params$sd^2,
    lognormal   = {
      s2 <- params$sdlog^2
      (exp(s2) - 1) * exp(2 * params$meanlog + s2)
    },
    uniform     = (params$max - params$min)^2 / 12,
    weibull     = {
      g1 <- gamma(1 + 1 / params$shape); g2 <- gamma(1 + 2 / params$shape)
      params$scale^2 * (g2 - g1^2)
    },
    tdist       = params$df / (params$df - 2),
    chisq       = 2 * params$df,
    fdist       = {
      d1 <- params$df1; d2 <- params$df2
      2 * d2^2 * (d1 + d2 - 2) / (d1 * (d2 - 2)^2 * (d2 - 4))
    },
    laplace     = 2 * params$scale^2,
    gumbel      = (pi^2 / 6) * params$scale^2,
    rayleigh    = params$scale^2 * (4 - pi) / 2,
    pareto      = {
      a <- params$shape; xm <- params$scale
      (xm^2 * a) / ((a - 1)^2 * (a - 2))
    },
    frechet     = {
      a <- params$shape; sg <- params$scale
      g1 <- gamma(1 - 1 / a); g2 <- gamma(1 - 2 / a)
      sg^2 * (g2 - g1^2)
    },
    inverse_gamma = {
      a <- params$shape; b <- params$scale
      b^2 / ((a - 1)^2 * (a - 2))
    },
    chi         = {
      m <- dist_mean_from_params("chi", params)
      params$df - m^2
    },
    folded_normal = {
      mu <- params$location; sigma <- params$scale
      m <- dist_mean_from_params("folded_normal", params)
      mu^2 + sigma^2 - m^2
    },
    erlang      = params$shape * params$scale^2,
    sym_triangular = params$scale^2 / 6,
    triangular  = {
      a <- params$a; b <- params$b; c <- params$c
      (a^2 + b^2 + c^2 - a * b - a * c - b * c) / 18
    },
    binomial    = params$size * params$prob * (1 - params$prob),
    poisson     = params$lambda,
    negative_binomial = params$size * (1 - params$prob) / params$prob^2,
    geometric   = (1 - params$prob) / params$prob^2,
    discrete_uniform = ((params$max - params$min + 1)^2 - 1) / 12,
    discrete_sym_triangular = params$n * (params$n + 2) / 6,
    discrete_triangular = {
      a <- params$a; b <- params$b; c <- params$c
      Z <- (b - a + 2) / 2
      ks <- a:b
      pmf <- ifelse(ks <= c, (ks - a + 1) / (c - a + 1) / Z,
                              (b - ks + 1) / (b - c + 1) / Z)
      m <- sum(ks * pmf)
      sum((ks - m)^2 * pmf)
    },
    stop(sprintf("no closed-form var wired for %s", name))
  )
}

# Canonical parameter list per family — used by partial_dist to identify the
# parameters a user may pin or leave free.
.CANONICAL_PARAMS <- list(
  gamma       = c("shape", "rate"),
  exponential = c("rate"),
  logistic    = c("location", "scale"),
  beta        = c("shape1", "shape2"),
  normal      = c("mean", "sd"),
  lognormal   = c("meanlog", "sdlog"),
  uniform     = c("min", "max"),
  weibull     = c("shape", "scale"),
  tdist       = c("df"),
  chisq       = c("df"),
  fdist       = c("df1", "df2"),
  laplace     = c("location", "scale"),
  gumbel      = c("location", "scale"),
  rayleigh    = c("scale"),
  pareto      = c("shape", "scale"),
  frechet     = c("shape", "scale"),
  inverse_gamma = c("shape", "scale"),
  chi         = c("df"),
  folded_normal = c("location", "scale"),
  erlang      = c("shape", "scale"),
  sym_triangular = c("location", "scale"),
  triangular  = c("a", "b", "c"),
  cauchy      = c("location", "scale"),
  binomial    = c("size", "prob"),
  poisson     = c("lambda"),
  negative_binomial = c("size", "prob"),
  geometric   = c("prob"),
  discrete_uniform = c("min", "max"),
  discrete_sym_triangular = c("mu", "n"),
  discrete_triangular = c("a", "b", "c")
)

canonical_params <- function(name) {
  p <- .CANONICAL_PARAMS[[name]]
  if (is.null(p)) stop(sprintf("No canonical params known for %s", name))
  p
}
