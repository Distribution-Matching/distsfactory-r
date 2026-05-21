# Discrete (asymmetric) triangular on integers {a, ..., b} with peak at c.
# PMF (Z = (b - a + 2)/2):
#   k in [a, c]:  P(k) = (k - a + 1) / (c - a + 1) / Z
#   k in [c+1, b]: P(k) = (b - k + 1) / (b - c + 1) / Z
# Identified by mean + var + mode: search around the continuous triangular
# solution within a ±1 integer neighborhood (mirrors the Python implementation).

.dtri_pmf_scalar <- function(k, a, b, c) {
  if (k < a || k > b) return(0)
  Z <- (b - a + 2) / 2
  if (k <= c) (k - a + 1) / (c - a + 1) / Z
  else        (b - k + 1) / (b - c + 1) / Z
}

ddtri_ <- function(x, a, b, c, log = FALSE) {
  out <- vapply(x, function(xi) {
    if (abs(xi - round(xi)) > 1e-9) return(0)
    .dtri_pmf_scalar(as.integer(round(xi)), a, b, c)
  }, numeric(1))
  if (log) log(out) else out
}

pdtri_ <- function(q, a, b, c, lower.tail = TRUE, log.p = FALSE) {
  out <- vapply(q, function(qi) {
    if (qi < a) return(0)
    if (qi >= b) return(1)
    sum(vapply(a:floor(qi), function(k) .dtri_pmf_scalar(k, a, b, c), numeric(1)))
  }, numeric(1))
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

qdtri_ <- function(p, a, b, c, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  vapply(p, function(pi_v) {
    if (pi_v <= 0) return(a)
    if (pi_v >= 1) return(b)
    s <- 0
    for (k in a:b) {
      s <- s + .dtri_pmf_scalar(k, a, b, c)
      if (s >= pi_v) return(k)
    }
    b
  }, numeric(1))
}

rdtri_ <- function(n, a, b, c) qdtri_(runif(n), a, b, c)

.dtri_mean <- function(a, b, c) {
  sum(vapply(a:b, function(k) k * .dtri_pmf_scalar(k, a, b, c), numeric(1)))
}

.dtri_var <- function(a, b, c) {
  m <- .dtri_mean(a, b, c)
  sum(vapply(a:b, function(k) (k - m)^2 * .dtri_pmf_scalar(k, a, b, c), numeric(1)))
}

dtri_from_mean_var_mode <- function(mean, var, mode) {
  if (var <= 0) stop("DiscreteTriangular variance must be positive")
  # Use the continuous triangular as a starting point.
  cont <- triang_from_mean_var_mode(mean, var, mode)
  a_cont <- cont$params$a; b_cont <- cont$params$b
  c_int <- as.integer(round(mode))
  a0 <- as.integer(round(a_cont))
  b0 <- as.integer(round(b_cont))

  best <- NULL
  best_err <- Inf
  for (da in c(-1, 0, 1)) {
    for (db in c(-1, 0, 1)) {
      a_try <- min(a0 + da, c_int)
      b_try <- max(b0 + db, c_int)
      if (!(a_try <= c_int && c_int <= b_try)) next
      m_try <- .dtri_mean(a_try, b_try, c_int)
      v_try <- .dtri_var(a_try, b_try, c_int)
      err <- ((m_try - mean) / max(abs(mean), 1))^2 +
             ((v_try - var)  / max(var, 1))^2
      if (err < best_err) {
        best_err <- err
        best <- list(a = a_try, b = b_try, c = c_int)
      }
    }
  }
  if (is.null(best))
    stop(sprintf("DiscreteTriangular: no integer triple near (mean=%g, var=%g, mode=%g)",
                 mean, var, mode))
  new_dist("discrete_triangular", best,
           function(x, a, b, c, log = FALSE) ddtri_(x, a, b, c, log),
           function(q, a, b, c, lower.tail = TRUE, log.p = FALSE)
             pdtri_(q, a, b, c, lower.tail, log.p),
           function(p, a, b, c, lower.tail = TRUE, log.p = FALSE)
             qdtri_(p, a, b, c, lower.tail, log.p),
           function(n, a, b, c) rdtri_(n, a, b, c))
}

dtri_exists_mean_var <- function(mean, var) {
  # Without a mode, infeasibility is undetermined; conservatively TRUE.
  is.finite(mean) && is.finite(var) && var > 0
}

dtri_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarModeSpec = dtri_from_mean_var_mode(spec$mean, spec$var, spec$mode),
    stop(sprintf("DiscreteTriangular does not support specification type '%s'", cls))
  )
}
