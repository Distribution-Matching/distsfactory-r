# F distribution matching (Fisher-Snedecor).
# Base R parameterization: df(x, df1, df2)
# mean = df2 / (df2 - 2)             for df2 > 2.
# var  = 2*df2^2 * (df1+df2-2) / (df1 * (df2-2)^2 * (df2-4))   for df2 > 4.
# Two free parameters (df1, df2); mean+var has a closed-form solution.

fdist_from_mean_var <- function(mean, var) {
  if (!is.finite(mean) || !is.finite(var)) stop("FDist mean/var must be finite")
  if (mean <= 1) stop(sprintf("FDist mean must exceed 1; got %g", mean))
  if (var <= 0) stop("FDist variance must be positive")
  # df2 from mean:  mean = df2/(df2-2)  ->  df2 = 2*mean/(mean-1)
  df2 <- 2 * mean / (mean - 1)
  if (df2 <= 4) {
    stop(sprintf("FDist mean=%g implies df2=%g <= 4 — variance undefined", mean, df2))
  }
  # var = 2*df2^2 * (df1+df2-2) / (df1 * (df2-2)^2 * (df2-4))
  # solve for df1: let C = 2*df2^2 / ((df2-2)^2 * (df2-4)).
  # var = C * (df1 + df2 - 2) / df1  ->  var * df1 = C*df1 + C*(df2-2)
  # df1 * (var - C) = C * (df2 - 2)
  # df1 = C * (df2 - 2) / (var - C)
  C <- 2 * df2^2 / ((df2 - 2)^2 * (df2 - 4))
  if (var <= C) {
    stop(sprintf("FDist mean=%g, var=%g: var must exceed %g", mean, var, C))
  }
  df1 <- C * (df2 - 2) / (var - C)
  if (df1 <= 0) stop("FDist mean+var: implied df1 must be positive")
  new_dist("fdist", list(df1 = df1, df2 = df2), df, pf, qf, rf)
}

fdist_exists_mean_var <- function(mean, var) {
  if (!is.finite(mean) || !is.finite(var) || mean <= 1 || var <= 0) return(FALSE)
  df2 <- 2 * mean / (mean - 1)
  if (df2 <= 4) return(FALSE)
  C <- 2 * df2^2 / ((df2 - 2)^2 * (df2 - 4))
  var > C
}

fdist_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = fdist_from_mean_var(spec$mean, spec$var),
    stop(sprintf("FDist does not support specification type '%s'", cls))
  )
}
