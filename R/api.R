# Public API: make_dist, dist_exists, available_distributions

# Registry of distribution handlers
.dist_handlers <- list(
  gamma = list(dispatch = gamma_dispatch, exists_mv = gamma_exists_mean_var),
  exponential = list(dispatch = exponential_dispatch, exists_mv = exponential_exists_mean_var),
  logistic = list(dispatch = logistic_dispatch, exists_mv = logistic_exists_mean_var),
  beta = list(dispatch = beta_dispatch, exists_mv = beta_exists_mean_var)
)

# Aliases
.dist_aliases <- c(exp = "exponential", expon = "exponential")

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
#' @param second_moment E[X^2] (requires \code{mean}).
#' @param median Target median (p = 0.5 quantile).
#' @param q1 First quartile (p = 0.25).
#' @param q3 Third quartile (p = 0.75).
#' @param iqr Interquartile range.
#' @param quantiles List of two \code{c(p, q)} vectors for arbitrary quantile
#'   constraints.
#' @param mode Target mode.
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
                      mode = NULL) {
  name <- resolve_dist_name(dist)
  handler <- .dist_handlers[[name]]
  spec <- parse_spec(
    mean = mean, var = var, std = std, cv = cv, scv = scv,
    second_moment = second_moment, median = median, q1 = q1, q3 = q3,
    iqr = iqr, quantiles = quantiles, mode = mode
  )
  handler$dispatch(spec)
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
                        mode = NULL) {
  name <- resolve_dist_name(dist)
  handler <- .dist_handlers[[name]]
  spec <- parse_spec(
    mean = mean, var = var, std = std, cv = cv, scv = scv,
    second_moment = second_moment, median = median, q1 = q1, q3 = q3,
    iqr = iqr, quantiles = quantiles, mode = mode
  )

  if (inherits(spec, "MeanVarSpec")) {
    return(handler$exists_mv(spec$mean, spec$var))
  }

  # For other specs, try construction
  tryCatch({
    make_dist(dist, mean = mean, var = var, std = std, cv = cv, scv = scv,
              second_moment = second_moment, median = median, q1 = q1,
              q3 = q3, iqr = iqr, quantiles = quantiles, mode = mode)
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
                                    quantiles = NULL, mode = NULL) {
  spec <- parse_spec(
    mean = mean, var = var, std = std, cv = cv, scv = scv,
    second_moment = second_moment, median = median, q1 = q1, q3 = q3,
    iqr = iqr, quantiles = quantiles, mode = mode
  )

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
