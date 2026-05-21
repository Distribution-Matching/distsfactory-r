# Symmetric Triangular distribution on [mu - a, mu + a].
# Parameters: location = mu, scale (half-width) = a > 0.
# pdf(x) = (a - |x - mu|) / a^2 on [mu-a, mu+a], else 0.
# mean   = mu;  var = a^2 / 6;  median = mode = mu.

dsymtri_ <- function(x, location = 0, scale = 1, log = FALSE) {
  if (scale <= 0) stop("dsymtri: scale must be positive")
  d <- abs(x - location)
  out <- ifelse(d >= scale, -Inf, log(scale - d) - 2 * log(scale))
  if (log) out else exp(out)
}

psymtri_ <- function(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("psymtri: scale must be positive")
  z <- (q - location) / scale  # in [-1, 1]
  p <- ifelse(z <= -1, 0,
       ifelse(z >=  1, 1,
       ifelse(z < 0, (1 + z)^2 / 2, 1 - (1 - z)^2 / 2)))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qsymtri_ <- function(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (scale <= 0) stop("qsymtri: scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  z <- ifelse(p < 0.5,
              sqrt(2 * p) - 1,
              1 - sqrt(2 * (1 - p)))
  location + scale * z
}

rsymtri_ <- function(n, location = 0, scale = 1) qsymtri_(runif(n), location, scale)

symtri_from_mean_var <- function(mean, var) {
  if (var <= 0) stop("SymTriangular variance must be positive")
  a <- sqrt(6 * var)
  new_dist("sym_triangular", list(location = mean, scale = a),
           dsymtri_, psymtri_, qsymtri_, rsymtri_)
}

symtri_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && var > 0
}

symtri_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec  = symtri_from_mean_var(spec$mean, spec$var),
    ModeVarSpec  = symtri_from_mean_var(spec$mode, spec$var),
    stop(sprintf("SymTriangular does not support specification type '%s'", cls))
  )
}
