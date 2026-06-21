# Distribution object: wraps base R d/p/q/r functions with baked-in parameters

#' Create a distribution object
#'
#' @param name Character name of the distribution (e.g. "gamma").
#' @param params Named list of parameters for the base R functions.
#' @param dfun,pfun,qfun,rfun The base R d/p/q/r functions.
#' @return A list with class "distsfactory_dist" containing d/p/q/r methods
#'   and a \code{$params} element.
#' @keywords internal
new_dist <- function(name, params, dfun, pfun, qfun, rfun) {
  obj <- list(
    name = name,
    params = params,
    d = function(x, ...) do.call(dfun, c(list(x), params, list(...))),
    p = function(q, ...) do.call(pfun, c(list(q), params, list(...))),
    q = function(p, ...) do.call(qfun, c(list(p), params, list(...))),
    r = function(n, ...) do.call(rfun, c(list(n), params, list(...))),
    mean = function() {
      tryCatch(dist_mean_from_params(name, params),
               error = function(e) NaN)
    },
    var = function() {
      tryCatch(dist_var_from_params(name, params),
               error = function(e) NaN)
    }
  )
  obj$std <- function() sqrt(obj$var())
  obj$median <- function() obj$q(0.5)
  class(obj) <- "distsfactory_dist"
  obj
}

#' Summary statistics for a distsfactory distribution
#'
#' S3 methods so that \code{mean(d)}, \code{median(d)}, and
#' \code{quantile(d, probs)} work on a \code{distsfactory_dist} just like
#' \code{d$mean()}, \code{d$median()}, \code{d$q(probs)}. Mirrors the Julia
#' (\code{mean(d)}, \code{var(d)}, \code{std(d)}, \code{median(d)}) and Python
#' (\code{d.mean()}, \code{d.var()}, \code{d.std()}, \code{d.median()})
#' sibling-package conventions. Use \code{d$var()} and \code{d$std()} for
#' variance and standard deviation since base R's \code{var()} is not generic.
#'
#' @param x A \code{distsfactory_dist} object.
#' @param probs Numeric vector of probabilities in \code{[0, 1]} (for
#'   \code{quantile}).
#' @param ... Unused.
#'
#' @return
#' \code{mean()} and \code{median()} each return a single numeric value of
#' length one (the mean and median of the distribution, respectively).
#' \code{quantile()} returns a named numeric vector of quantiles, one element
#' per probability in \code{probs}, with names giving the corresponding
#' percentages.
#'
#' @export
mean.distsfactory_dist <- function(x, ...) x$mean()

#' @rdname mean.distsfactory_dist
#' @export
median.distsfactory_dist <- function(x, ...) x$median()

#' @rdname mean.distsfactory_dist
#' @export
quantile.distsfactory_dist <- function(x, probs = seq(0, 1, 0.25), ...) {
  out <- x$q(probs)
  names(out) <- paste0(round(probs * 100, 4), "%")
  out
}

#' Print a distsfactory distribution
#'
#' @param x A \code{distsfactory_dist} object.
#' @param ... Unused.
#'
#' @return Invisibly returns the \code{distsfactory_dist} object \code{x}.
#'   Called for its side effect of printing a description of the distribution
#'   to the console.
#'
#' @export
print.distsfactory_dist <- function(x, ...) {
  fmt <- function(v) {
    if (length(v) > 1) {
      sprintf("c(%s)", paste(format(v, digits = 4), collapse = ", "))
    } else {
      format(v, digits = 4)
    }
  }
  if (!is.null(x$params)) {
    # Plain (un-wrapped) family dist: show its canonical parameters.
    param_str <- paste(names(x$params), "=",
                       vapply(x$params, fmt, character(1)),
                       collapse = ", ")
    cat(sprintf("distsfactory: %s(%s)\n", x$name, param_str))
  } else {
    # Wrapped dist (truncated / shifted / flipped / scaled): describe the
    # transform and the parent. No $params slot (the parent's parameters are
    # distinct from the wrapped dist's moments).
    sup <- x$support
    sup_str <- if (!is.null(sup)) sprintf("[%g, %g]", sup[1], sup[2]) else "(no support)"
    transform <- if (!is.null(x$shift))      sprintf("shifted by %g",    x$shift)
                 else if (!is.null(x$flip_point)) sprintf("flipped about %g", x$flip_point)
                 else if (!is.null(x$scale_loc))  sprintf("scaled to %s",     sup_str)
                 else                                sprintf("truncated to %s",  sup_str)
    parent_str <- paste(names(x$parent$params), "=",
                        vapply(x$parent$params, fmt, character(1)),
                        collapse = ", ")
    cat(sprintf("distsfactory: %s %s\n  parent: %s(%s)\n",
                x$name, transform, x$parent$name, parent_str))
  }
  invisible(x)
}
