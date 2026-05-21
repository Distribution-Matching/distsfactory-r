# Discrete symmetric triangular on {mu - n, ..., mu + n}.
# PMF: P(mu + k) = (n + 1 - |k|) / (n + 1)^2  for k in {-n, ..., n}
# mean = mu;  var = n*(n+2)/6;  mode = mu.

.dst_size <- function(n) (n + 1)^2

ddst_ <- function(x, mu = 0, n = 1, log = FALSE) {
  if (n < 0 || abs(n - round(n)) > 1e-9) stop("ddst_: n must be a non-negative integer")
  k <- x - mu
  out <- ifelse(abs(k) > n | abs(k - round(k)) > 1e-9,
                -Inf,
                log(n + 1 - abs(k)) - 2 * log(n + 1))
  if (log) out else exp(out)
}

pdst_ <- function(q, mu = 0, n = 1, lower.tail = TRUE, log.p = FALSE) {
  if (n < 0) stop("pdst_: n must be >= 0")
  Z <- (n + 1)^2
  out <- numeric(length(q))
  for (i in seq_along(q)) {
    qi <- q[i]
    if (qi < mu - n) { out[i] <- 0; next }
    if (qi >= mu + n) { out[i] <- 1; next }
    k_max <- floor(qi) - mu
    s <- 0
    for (k in seq.int(-n, k_max)) s <- s + (n + 1 - abs(k))
    out[i] <- s / Z
  }
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

qdst_ <- function(p, mu = 0, n = 1, lower.tail = TRUE, log.p = FALSE) {
  if (n < 0) stop("qdst_: n must be >= 0")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  Z <- (n + 1)^2
  out <- numeric(length(p))
  for (i in seq_along(p)) {
    pi_v <- p[i]
    if (pi_v <= 0) { out[i] <- mu - n; next }
    if (pi_v >= 1) { out[i] <- mu + n; next }
    s <- 0; chosen <- mu + n
    for (k in seq.int(-n, n)) {
      s <- s + (n + 1 - abs(k)) / Z
      if (s >= pi_v) { chosen <- mu + k; break }
    }
    out[i] <- chosen
  }
  out
}

rdst_ <- function(N, mu = 0, n = 1) qdst_(runif(N), mu = mu, n = n)

dst_from_mean_var <- function(mean, var) {
  if (var <= 0) stop("DiscreteSymTriangular variance must be positive")
  # var = n*(n+2)/6  ->  n = -1 + sqrt(1 + 6*var)
  n_real <- -1 + sqrt(1 + 6 * var)
  n <- round(n_real)
  if (n < 1) stop("DiscreteSymTriangular: implied n < 1")
  mu <- round(mean)
  new_dist("discrete_sym_triangular", list(mu = mu, n = n),
           ddst_, pdst_, qdst_, rdst_)
}

dst_exists_mean_var <- function(mean, var) {
  is.finite(mean) && is.finite(var) && var > 0 &&
    round(-1 + sqrt(1 + 6 * var)) >= 1
}

dst_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarSpec = dst_from_mean_var(spec$mean, spec$var),
    stop(sprintf("DiscreteSymTriangular does not support specification type '%s'", cls))
  )
}
