# (Asymmetric) Triangular distribution on [a, b] with peak at c, a <= c <= b.
# Three parameters; identified by mean + var + mode.
# mean = (a + b + c) / 3
# var  = (a^2 + b^2 + c^2 - a*b - a*c - b*c) / 18
# mode = c
# pdf(x) =
#   2(x-a) / ((b-a)(c-a))   for a <= x <= c
#   2(b-x) / ((b-a)(b-c))   for c <= x <= b
# cdf(x) =
#   (x-a)^2 / ((b-a)(c-a))  for a <= x <= c
#   1 - (b-x)^2 / ((b-a)(b-c)) for c < x <= b

dtriang_ <- function(x, a, b, c, log = FALSE) {
  if (!(a < b) || !(a <= c && c <= b)) stop("dtriang: require a < b and a <= c <= b")
  out <- ifelse(x < a | x > b, -Inf,
         ifelse(x <= c,
                log(2) + log(x - a) - log(b - a) - log(c - a),
                log(2) + log(b - x) - log(b - a) - log(b - c)))
  if (log) out else exp(out)
}

ptriang_ <- function(q, a, b, c, lower.tail = TRUE, log.p = FALSE) {
  if (!(a < b) || !(a <= c && c <= b)) stop("ptriang: invalid (a,b,c)")
  p <- ifelse(q <= a, 0,
       ifelse(q >= b, 1,
       ifelse(q <= c, (q - a)^2 / ((b - a) * (c - a)),
                       1 - (b - q)^2 / ((b - a) * (b - c)))))
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

qtriang_ <- function(p, a, b, c, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  cp <- (c - a) / (b - a)
  ifelse(p <= cp,
         a + sqrt(p * (b - a) * (c - a)),
         b - sqrt((1 - p) * (b - a) * (b - c)))
}

rtriang_ <- function(n, a, b, c) qtriang_(runif(n), a, b, c)

triang_from_mean_var_mode <- function(mean, var, mode) {
  # mean = (a+b+c)/3 with c = mode -> a + b = 3*mean - mode  =: S
  # var formula in (a, b) given c = mode:
  #   18*var = a^2 + b^2 + c^2 - ab - ac - bc
  # Let S = a+b, P = a*b. Then a^2+b^2 = S^2 - 2P; ab = P.
  # 18*var = S^2 - 2P + c^2 - P - c*(a+b) = S^2 - 3P + c^2 - c*S
  # => P = (S^2 - c*S + c^2 - 18*var) / 3
  S <- 3 * mean - mode
  c0 <- mode
  P <- (S^2 - c0 * S + c0^2 - 18 * var) / 3
  # a, b are roots of t^2 - S*t + P = 0
  disc <- S^2 - 4 * P
  if (disc < 0) stop("Triangular mean+var+mode infeasible (negative discriminant)")
  a <- (S - sqrt(disc)) / 2
  b <- (S + sqrt(disc)) / 2
  if (!(a <= c0 && c0 <= b))
    stop(sprintf("Triangular: implied (a=%g, b=%g, c=%g) violates a <= c <= b", a, b, c0))
  new_dist("triangular", list(a = a, b = b, c = c0),
           function(x, a, b, c, log = FALSE) dtriang_(x, a, b, c, log),
           function(q, a, b, c, lower.tail = TRUE, log.p = FALSE)
             ptriang_(q, a, b, c, lower.tail, log.p),
           function(p, a, b, c, lower.tail = TRUE, log.p = FALSE)
             qtriang_(p, a, b, c, lower.tail, log.p),
           function(n, a, b, c) rtriang_(n, a, b, c))
}

triang_exists_mean_var <- function(mean, var) {
  # Without a mode, infeasibility is undetermined; conservatively say TRUE if
  # var > 0 and mean is finite. This is consistent with available_distributions
  # offering it when given enough specs.
  is.finite(mean) && is.finite(var) && var > 0
}

# Triangular dispatch is unusual: it takes mode+var+mean, which parse_spec
# does not produce on its own. Add an explicit hook by extending parse_spec
# below — see R/spec.R for the MeanVarModeSpec class.
triang_dispatch <- function(spec) {
  cls <- class(spec)
  switch(cls,
    MeanVarModeSpec = triang_from_mean_var_mode(spec$mean, spec$var, spec$mode),
    stop(sprintf("Triangular does not support specification type '%s'", cls))
  )
}
