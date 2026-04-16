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
    r = function(n, ...) do.call(rfun, c(list(n), params, list(...)))
  )
  class(obj) <- "distsfactory_dist"
  obj
}

#' @export
print.distsfactory_dist <- function(x, ...) {
  param_str <- paste(names(x$params), "=", format(unlist(x$params), digits = 4),
                     collapse = ", ")
  cat(sprintf("distsfactory: %s(%s)\n", x$name, param_str))
  invisible(x)
}
