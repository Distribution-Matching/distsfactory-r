# Discrete uniform on integers {a, a+1, ..., b}.
# mean = (a+b)/2
# var  = ((b-a+1)^2 - 1) / 12
# Equivalently with n = b - a (so |support| = n+1):  var = n*(n+2)/12.

ddunif_ <- function(x, min = 0, max = 1, log = FALSE) {
  size <- max - min + 1
  if (size <= 0) stop("ddunif: max must be >= min")
  out <- ifelse(x < min | x > max | abs(x - round(x)) > 1e-9, -Inf, -log(size))
  if (log) out else exp(out)
}

pdunif_ <- function(q, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE) {
  size <- max - min + 1
  if (size <= 0) stop("pdunif: max must be >= min")
  k <- floor(q)
  p <- ifelse(q < min, 0,
       ifelse(q >= max, 1, (k - min + 1) / size))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qdunif_ <- function(p, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE) {
  size <- max - min + 1
  if (size <= 0) stop("qdunif: max must be >= min")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  k <- ceiling(p * size) + min - 1
  pmin(pmax(k, min), max)
}

rdunif_ <- function(n, min = 0, max = 1) {
  sample.int(max - min + 1, size = n, replace = TRUE) + min - 1
}

dunif_from_mean_var <- function(mean, var) {
  if (var <= 0) stop("DiscreteUniform variance must be positive")
  # var = ((b - a + 1)^2 - 1) / 12  ->  b - a + 1 = sqrt(12*var + 1)
  span <- sqrt(12 * var + 1)
  n <- round(span - 1)  # b - a = n
  if (n < 0) stop("DiscreteUniform: implied support empty")
  a <- round(mean - n / 2)
  b <- a + n
  new_dist("discrete_uniform", list(min = a, max = b),
           ddunif_, pdunif_, qdunif_, rdunif_)
}

dunif_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && var > 0
}

dunif_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = dunif_from_mean_var(spec$mean, spec$var),
    stop(sprintf("DiscreteUniform does not support specification type '%s'", cls))
  )
}
