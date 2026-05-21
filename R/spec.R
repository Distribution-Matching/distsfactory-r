# Specification parser: converts keyword arguments into typed spec lists

#' Parse keyword arguments into a specification object
#' @param mean Target mean.
#' @param var Target variance.
#' @param std Target standard deviation.
#' @param cv Coefficient of variation.
#' @param scv Squared coefficient of variation.
#' @param second_moment Target second raw moment.
#' @param median Target median.
#' @param q1 First quartile.
#' @param q3 Third quartile.
#' @param iqr Interquartile range.
#' @param quantiles List of two \code{c(p, q)} vectors.
#' @param mode Target mode.
#' @return A list with class set to the spec type.
#' @keywords internal
parse_spec <- function(mean = NULL, var = NULL, std = NULL, cv = NULL,
                       scv = NULL, second_moment = NULL, median = NULL,
                       q1 = NULL, q3 = NULL, iqr = NULL, quantiles = NULL,
                       mode = NULL) {
  # Resolve variance from alternative dispersion measures
  if (is.null(var) && !is.null(std)) {
    var <- std^2
  } else if (is.null(var) && !is.null(mean) && !is.null(cv)) {
    var <- (cv * mean)^2
  } else if (is.null(var) && !is.null(mean) && !is.null(scv)) {
    var <- scv * mean^2
  } else if (is.null(var) && !is.null(mean) && !is.null(second_moment)) {
    var <- second_moment - mean^2
  }

  # Mode-based specs
  if (!is.null(mode)) {
    if (!is.null(mean) && !is.null(var))
      return(structure(list(mean = mean, var = var, mode = mode),
                       class = "MeanVarModeSpec"))
    if (!is.null(mean) && is.null(var))
      return(structure(list(mean = mean, mode = mode), class = "MeanModeSpec"))
    if (!is.null(var))
      return(structure(list(mode = mode, var = var), class = "ModeVarSpec"))
    if (!is.null(median))
      return(structure(list(mode = mode, p = 0.5, q = median), class = "ModeQuantileSpec"))
    if (!is.null(q1))
      return(structure(list(mode = mode, p = 0.25, q = q1), class = "ModeQuantileSpec"))
    if (!is.null(q3))
      return(structure(list(mode = mode, p = 0.75, q = q3), class = "ModeQuantileSpec"))
    if (!is.null(iqr))
      return(structure(list(mode = mode, iqr = iqr), class = "ModeIQRSpec"))
    # Mode alone — only some 1-parameter families accept this (e.g. Rayleigh).
    return(structure(list(mode = mode), class = "ModeSpec"))
  }

  # Quantile-based specs
  if (!is.null(quantiles)) {
    if (length(quantiles) != 2) stop("quantiles must be a list of 2 c(p, q) vectors")
    return(structure(list(
      p1 = quantiles[[1]][1], q1 = quantiles[[1]][2],
      p2 = quantiles[[2]][1], q2 = quantiles[[2]][2]
    ), class = "TwoQuantileSpec"))
  }
  if (!is.null(q1) && !is.null(q3))
    return(structure(list(p1 = 0.25, q1 = q1, p2 = 0.75, q2 = q3), class = "TwoQuantileSpec"))
  if (!is.null(median) && !is.null(iqr))
    return(structure(list(
      p1 = 0.25, q1 = median - iqr / 2,
      p2 = 0.75, q2 = median + iqr / 2
    ), class = "TwoQuantileSpec"))
  if (!is.null(mean) && !is.null(median))
    return(structure(list(mean = mean, p = 0.5, q = median), class = "MeanQuantileSpec"))
  if (!is.null(mean) && !is.null(q1))
    return(structure(list(mean = mean, p = 0.25, q = q1), class = "MeanQuantileSpec"))
  if (!is.null(mean) && !is.null(q3))
    return(structure(list(mean = mean, p = 0.75, q = q3), class = "MeanQuantileSpec"))
  if (!is.null(median))
    return(structure(list(p = 0.5, q = median), class = "QuantileSpec"))
  if (!is.null(q1))
    return(structure(list(p = 0.25, q = q1), class = "QuantileSpec"))
  if (!is.null(q3))
    return(structure(list(p = 0.75, q = q3), class = "QuantileSpec"))

  # Moment-based specs
  if (!is.null(mean) && !is.null(var))
    return(structure(list(mean = mean, var = var), class = "MeanVarSpec"))
  if (!is.null(mean))
    return(structure(list(mean = mean), class = "MeanSpec"))
  if (!is.null(var))
    return(structure(list(var = var), class = "VarSpec"))

  stop("Must provide at least one moment or quantile specification")
}
