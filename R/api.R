# Public API: make_dist, dist_exists, available_distributions

# Registry of distribution handlers
.dist_handlers <- list(
  gamma       = list(dispatch = gamma_dispatch,       exists_mv = gamma_exists_mean_var),
  exponential = list(dispatch = exponential_dispatch, exists_mv = exponential_exists_mean_var),
  logistic    = list(dispatch = logistic_dispatch,    exists_mv = logistic_exists_mean_var),
  beta        = list(dispatch = beta_dispatch,        exists_mv = beta_exists_mean_var),
  normal      = list(dispatch = normal_dispatch,      exists_mv = normal_exists_mean_var),
  lognormal   = list(dispatch = lognormal_dispatch,   exists_mv = lognormal_exists_mean_var),
  uniform     = list(dispatch = uniform_dispatch,     exists_mv = uniform_exists_mean_var),
  weibull     = list(dispatch = weibull_dispatch,     exists_mv = weibull_exists_mean_var),
  tdist       = list(dispatch = tdist_dispatch,       exists_mv = tdist_exists_mean_var),
  chisq       = list(dispatch = chisq_dispatch,       exists_mv = chisq_exists_mean_var),
  fdist       = list(dispatch = fdist_dispatch,       exists_mv = fdist_exists_mean_var),
  laplace     = list(dispatch = laplace_dispatch,     exists_mv = laplace_exists_mean_var),
  gumbel      = list(dispatch = gumbel_dispatch,      exists_mv = gumbel_exists_mean_var),
  rayleigh    = list(dispatch = rayleigh_dispatch,    exists_mv = rayleigh_exists_mean_var),
  pareto      = list(dispatch = pareto_dispatch,      exists_mv = pareto_exists_mean_var),
  frechet     = list(dispatch = frechet_dispatch,     exists_mv = frechet_exists_mean_var),
  inverse_gamma = list(dispatch = invgamma_dispatch,  exists_mv = invgamma_exists_mean_var),
  chi         = list(dispatch = chi_dispatch,         exists_mv = chi_exists_mean_var),
  folded_normal = list(dispatch = foldednorm_dispatch, exists_mv = foldednorm_exists_mean_var),
  erlang      = list(dispatch = erlang_dispatch,      exists_mv = erlang_exists_mean_var),
  sym_triangular = list(dispatch = symtri_dispatch,   exists_mv = symtri_exists_mean_var),
  triangular  = list(dispatch = triang_dispatch,      exists_mv = triang_exists_mean_var),
  cauchy      = list(dispatch = cauchy_dispatch,      exists_mv = cauchy_exists_mean_var),
  binomial    = list(dispatch = binomial_dispatch,    exists_mv = binomial_exists_mean_var),
  poisson     = list(dispatch = poisson_dispatch,     exists_mv = poisson_exists_mean_var),
  negative_binomial = list(dispatch = nbinom_dispatch, exists_mv = nbinom_exists_mean_var),
  geometric   = list(dispatch = geometric_dispatch,   exists_mv = geometric_exists_mean_var),
  discrete_uniform = list(dispatch = dunif_dispatch,  exists_mv = dunif_exists_mean_var),
  discrete_sym_triangular = list(dispatch = dst_dispatch, exists_mv = dst_exists_mean_var),
  discrete_triangular = list(dispatch = dtri_dispatch, exists_mv = dtri_exists_mean_var)
)

# Aliases
.dist_aliases <- c(
  exp       = "exponential",
  expon     = "exponential",
  norm      = "normal",
  gauss     = "normal",
  gaussian  = "normal",
  lnorm     = "lognormal",
  log_normal = "lognormal",
  unif      = "uniform",
  weib      = "weibull",
  student   = "tdist",
  t         = "tdist",
  chi_sq    = "chisq",
  chisquare = "chisq",
  chi_squared = "chisq",
  f         = "fdist",
  fisher    = "fdist",
  invgamma  = "inverse_gamma",
  "inverse-gamma" = "inverse_gamma",
  foldednormal = "folded_normal",
  symtriangular = "sym_triangular",
  triang = "triangular",
  nbinom = "negative_binomial",
  "negbinom" = "negative_binomial",
  binom = "binomial",
  geom = "geometric",
  pois = "poisson",
  dunif = "discrete_uniform",
  randint = "discrete_uniform"
)

resolve_dist_name <- function(dist) {
  name <- tolower(dist)
  if (name %in% names(.dist_aliases)) name <- .dist_aliases[[name]]
  if (!name %in% names(.dist_handlers)) {
    stop(sprintf(
      "Unknown distribution: '%s'. Available: %s",
      dist, paste(sort(names(.dist_handlers)), collapse = ", ")
    ))
  }
  name
}


#' Construct a distribution from partial specifications
#'
#' Creates a distribution object whose moments or quantiles match the
#' given constraints. Returns an object with \code{$d()}, \code{$p()},
#' \code{$q()}, \code{$r()} methods wrapping base R distribution functions,
#' plus a \code{$params} list.
#'
#' @param dist Character string naming the distribution (case-insensitive).
#'   Supported: \code{"gamma"}, \code{"exponential"} (or \code{"exp"}),
#'   \code{"logistic"}, \code{"beta"}.
#' @param mean Target mean.
#' @param var Target variance.
#' @param std Target standard deviation (converted to variance).
#' @param cv Coefficient of variation (requires \code{mean}).
#' @param scv Squared coefficient of variation (requires \code{mean}).
#' @param second_moment Target second raw moment (requires \code{mean}).
#' @param median Target median (p = 0.5 quantile).
#' @param q1 First quartile (p = 0.25).
#' @param q3 Third quartile (p = 0.75).
#' @param iqr Interquartile range.
#' @param quantiles List of two \code{c(p, q)} vectors for arbitrary quantile
#'   constraints.
#' @param mode Target mode.
#' @param support Optional length-2 numeric vector \code{c(lo, hi)} giving the
#'   target support; either endpoint may be \code{Inf}. When given, the
#'   distribution is placed on this support via an affine transform (when the
#'   requested shape matches the natural one) or truncation (when it is
#'   strictly smaller). Currently used together with \code{mean} and
#'   \code{var}.
#'
#' @return An object of class \code{"distsfactory_dist"} with elements:
#'   \describe{
#'     \item{\code{$d(x)}}{Density function}
#'     \item{\code{$p(q)}}{CDF}
#'     \item{\code{$q(p)}}{Quantile function}
#'     \item{\code{$r(n)}}{Random generation}
#'     \item{\code{$params}}{Named list of distribution parameters}
#'     \item{\code{$name}}{Distribution name}
#'   }
#'
#' @examples
#' d <- make_dist("gamma", mean = 5, var = 3)
#' d$d(2)       # density at x = 2
#' d$p(5)       # CDF at x = 5
#' d$q(0.5)     # median
#' d$r(10)      # 10 random samples
#' d$params     # list(shape = ..., rate = ...)
#'
#' @export
make_dist <- function(dist, mean = NULL, var = NULL, std = NULL, cv = NULL,
                      scv = NULL, second_moment = NULL, median = NULL,
                      q1 = NULL, q3 = NULL, iqr = NULL, quantiles = NULL,
                      mode = NULL, support = NULL) {
  if (inherits(dist, "partial_dist")) {
    return(make_dist_from_partial(dist, mean = mean, var = var, std = std,
                                  cv = cv, scv = scv,
                                  second_moment = second_moment))
  }
  name <- resolve_dist_name(dist)
  handler <- .dist_handlers[[name]]
  spec <- parse_spec(
    mean = mean, var = var, std = std, cv = cv, scv = scv,
    second_moment = second_moment, median = median, q1 = q1, q3 = q3,
    iqr = iqr, quantiles = quantiles, mode = mode
  )
  if (is.null(support)) {
    return(handler$dispatch(spec))
  }
  dist_on_support(name, spec, support)
}


#' Check if a distribution can be constructed with given constraints
#'
#' @inheritParams make_dist
#' @return \code{TRUE} if feasible, \code{FALSE} otherwise.
#'
#' @examples
#' dist_exists("beta", mean = 0.5, var = 0.1)        # TRUE
#' dist_exists("beta", mean = 0.5, var = 0.3)        # FALSE
#' dist_exists("exponential", mean = 2.5, var = 6.25) # TRUE
#'
#' @export
dist_exists <- function(dist, mean = NULL, var = NULL, std = NULL, cv = NULL,
                        scv = NULL, second_moment = NULL, median = NULL,
                        q1 = NULL, q3 = NULL, iqr = NULL, quantiles = NULL,
                        mode = NULL, support = NULL) {
  name <- resolve_dist_name(dist)
  handler <- .dist_handlers[[name]]
  spec <- parse_spec(
    mean = mean, var = var, std = std, cv = cv, scv = scv,
    second_moment = second_moment, median = median, q1 = q1, q3 = q3,
    iqr = iqr, quantiles = quantiles, mode = mode
  )

  if (!is.null(support) && inherits(spec, "MeanVarSpec")) {
    return(dist_exists_on_support(name, spec$mean, spec$var, support))
  }

  if (inherits(spec, "MeanVarSpec")) {
    return(handler$exists_mv(spec$mean, spec$var))
  }

  # For other specs, try construction
  tryCatch({
    make_dist(dist, mean = mean, var = var, std = std, cv = cv, scv = scv,
              second_moment = second_moment, median = median, q1 = q1,
              q3 = q3, iqr = iqr, quantiles = quantiles, mode = mode,
              support = support)
    TRUE
  }, error = function(e) FALSE)
}


#' List distributions feasible for the given constraints
#'
#' @inheritParams make_dist
#' @return Character vector of distribution names.
#'
#' @examples
#' available_distributions(mean = 5, var = 3)
#' available_distributions(mean = 0.5, var = 0.05)
#'
#' @export
available_distributions <- function(mean = NULL, var = NULL, std = NULL,
                                    cv = NULL, scv = NULL,
                                    second_moment = NULL, median = NULL,
                                    q1 = NULL, q3 = NULL, iqr = NULL,
                                    quantiles = NULL, mode = NULL,
                                    support = NULL) {
  if (all(vapply(list(mean, var, std, cv, scv, second_moment, median,
                      q1, q3, iqr, quantiles, mode, support),
                 is.null, logical(1)))) {
    # No constraints — return the full registered set.
    return(names(.dist_handlers))
  }
  spec <- parse_spec(
    mean = mean, var = var, std = std, cv = cv, scv = scv,
    second_moment = second_moment, median = median, q1 = q1, q3 = q3,
    iqr = iqr, quantiles = quantiles, mode = mode
  )

  if (!is.null(support) && inherits(spec, "MeanVarSpec")) {
    feasible <- character(0)
    for (name in names(.dist_handlers)) {
      ok <- tryCatch(
        dist_exists_on_support(name, spec$mean, spec$var, support),
        error = function(e) FALSE
      )
      if (isTRUE(ok)) feasible <- c(feasible, name)
    }
    return(feasible)
  }

  if (inherits(spec, "MeanVarSpec")) {
    feasible <- character(0)
    for (name in names(.dist_handlers)) {
      if (.dist_handlers[[name]]$exists_mv(spec$mean, spec$var)) {
        feasible <- c(feasible, name)
      }
    }
    return(feasible)
  }

  # For other specs, try construction
  feasible <- character(0)
  for (name in names(.dist_handlers)) {
    ok <- tryCatch({
      make_dist(name, mean = mean, var = var, std = std, cv = cv, scv = scv,
                second_moment = second_moment, median = median, q1 = q1,
                q3 = q3, iqr = iqr, quantiles = quantiles, mode = mode)
      TRUE
    }, error = function(e) FALSE)
    if (ok) feasible <- c(feasible, name)
  }
  feasible
}
