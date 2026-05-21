# Construct / check distributions on arbitrary supports.
#
# Mirrors src/support.jl in DistributionsFactories.jl. Given a target support,
# chooses between:
#
# - **Affine transform** when the requested support has the same shape as the
#   distribution's natural support (e.g. Gamma on [a, Inf), Beta on [a, b]).
# - **Truncation** when the requested support is strictly contained in the
#   natural one (e.g. Normal on [a, b]).
#
# Continuous supports are passed as a 2-vector c(lo, hi); either endpoint may
# be +Inf or -Inf. Discrete supports are passed as an integer 2-vector
# c(a, b) representing the inclusive integer interval {a, ..., b}.

# Natural support classification by canonical family name.
.SUPPORT_TYPE <- c(
  normal            = "real",
  laplace           = "real",
  logistic          = "real",
  gumbel            = "real",
  tdist             = "real",
  uniform           = "real",
  sym_triangular    = "real",
  triangular        = "real",
  gamma             = "positive",
  erlang            = "positive",
  exponential       = "positive",
  lognormal         = "positive",
  weibull           = "positive",
  frechet           = "positive",
  chi               = "positive",
  chisq             = "positive",
  rayleigh          = "positive",
  fdist             = "positive",
  inverse_gamma     = "positive",
  pareto            = "positive",
  folded_normal     = "positive",
  beta              = "unit",
  binomial          = "integer_bounded",
  discrete_uniform  = "integer_bounded",
  discrete_sym_triangular = "integer_bounded",
  discrete_triangular = "integer_bounded",
  poisson           = "integer_nonneg",
  negative_binomial = "integer_nonneg",
  geometric         = "integer_nonneg"
)

support_endpoints <- function(support) {
  if (length(support) != 2) stop("support must be a length-2 vector c(lo, hi)")
  lo <- support[[1]]; hi <- support[[2]]
  list(lo = as.numeric(lo), hi = as.numeric(hi))
}

requested_support_shape <- function(lo, hi) {
  if (is.infinite(lo) && lo < 0 && is.infinite(hi) && hi > 0) return("real")
  if (is.finite(lo) && is.infinite(hi) && hi > 0) return("half_right")
  if (is.infinite(lo) && lo < 0 && is.finite(hi)) return("half_left")
  if (is.finite(lo) && is.finite(hi)) return("bounded")
  stop(sprintf("Invalid support endpoints: (%g, %g)", lo, hi))
}

natural_support_of <- function(name) {
  s <- .SUPPORT_TYPE[[name]]
  if (is.null(s)) stop(sprintf("%s: no natural support classification known", name))
  s
}


# ---------------------------------------------------------------------------
# Feasibility predicate on arbitrary supports.
# ---------------------------------------------------------------------------

dist_exists_on_support <- function(name, mean, var, support) {
  # Mirrors Julia's `_dist_exists_on_support` and Python's
  # `why_not_mean_var_on_support`: the Langevin envelope is *not* applied here
  # for truncated locscale families — only the structural natural-family check
  # after standardizing moments to the natural support. The tighter Langevin
  # dome is enforced inside the constructor (`make_dist(..., support=)`).
  ep <- support_endpoints(support)
  lo <- ep$lo; hi <- ep$hi
  natural <- natural_support_of(name)
  requested <- requested_support_shape(lo, hi)
  exists_mv <- .dist_handlers[[name]]$exists_mv

  if (natural == "real") {
    return(exists_mv(mean, var))
  }

  if (natural == "positive") {
    if (requested == "half_right") return(exists_mv(mean - lo, var))
    if (requested == "half_left")  return(exists_mv(hi - mean, var))
    if (requested == "bounded") {
      if (lo < 0) return(FALSE)
      return(exists_mv(mean, var))
    }
    return(FALSE)  # natural=positive cannot be placed on (-Inf, Inf)
  }

  if (natural == "unit") {
    if (requested == "bounded") {
      w <- hi - lo
      return(exists_mv((mean - lo) / w, var / w^2))
    }
    return(FALSE)
  }

  # Discrete supports — defer to natural predicate.
  exists_mv(mean, var)
}


# ---------------------------------------------------------------------------
# Truncation wrapper: build a frozen-like dist whose d/p/q/r are the inner
# dist conditioned on x in [lo, hi].
# ---------------------------------------------------------------------------

make_truncated <- function(inner_name, inner_params, dfun, pfun, qfun, rfun, lo, hi) {
  F_lo <- if (is.infinite(lo) && lo < 0) 0
          else do.call(pfun, c(list(lo), inner_params))
  F_hi <- if (is.infinite(hi) && hi > 0) 1
          else do.call(pfun, c(list(hi), inner_params))
  Z <- F_hi - F_lo
  if (!(Z > 0))
    stop(sprintf("Truncated %s: truncation mass is zero on [%g, %g]", inner_name, lo, hi))

  params <- c(inner_params, list(.lo = lo, .hi = hi))

  d_trunc <- function(x, log = FALSE) {
    raw <- do.call(dfun, c(list(x), inner_params))
    in_support <- !(x < lo | x > hi)
    if (log) {
      log_raw <- do.call(dfun, c(list(x), inner_params, list(log = TRUE)))
      ifelse(in_support, log_raw - log(Z), -Inf)
    } else {
      ifelse(in_support, raw / Z, 0)
    }
  }
  p_trunc <- function(q, lower.tail = TRUE, log.p = FALSE) {
    raw <- do.call(pfun, c(list(pmin(pmax(q, lo), hi)), inner_params))
    p <- (raw - F_lo) / Z
    if (!lower.tail) p <- 1 - p
    if (log.p) log(p) else p
  }
  q_trunc <- function(p, lower.tail = TRUE, log.p = FALSE) {
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    do.call(qfun, c(list(F_lo + p * Z), inner_params))
  }
  r_trunc <- function(n) {
    u <- runif(n)
    do.call(qfun, c(list(F_lo + u * Z), inner_params))
  }

  out <- list(
    name = paste0("truncated_", inner_name),
    params = params,
    d = d_trunc, p = p_trunc, q = q_trunc, r = r_trunc
  )
  class(out) <- "distsfactory_dist"
  out
}


# ---------------------------------------------------------------------------
# Constructor on arbitrary support — make_dist(... support=).
# ---------------------------------------------------------------------------

dist_on_support <- function(name, spec, support) {
  ep <- support_endpoints(support)
  lo <- ep$lo; hi <- ep$hi
  natural <- natural_support_of(name)
  requested <- requested_support_shape(lo, hi)

  if (requested == "real" && natural == "real") {
    return(.dist_handlers[[name]]$dispatch(spec))
  }

  if (inherits(spec, "MeanVarSpec")) {
    mu <- spec$mean; var <- spec$var
  } else {
    stop("support= currently only supports MeanVarSpec (mean+var) constructions")
  }

  if (natural == "real") {
    if (name %in% c("normal", "laplace", "logistic")) {
      if (!truncated_locscale_exists(name, mu, var, lo, hi))
        stop(sprintf("Truncated %s on [%g, %g] infeasible at (mean=%g, var=%g)",
                     name, lo, hi, mu, var))
      return(solve_truncated_locscale(name, mu, var, lo, hi))
    }
    stop(sprintf("Truncation of %s on [%g, %g] not implemented", name, lo, hi))
  }

  if (natural == "positive") {
    if (requested == "half_right") {
      # Affine shift onto [lo, Inf)
      inner <- .dist_handlers[[name]]$dispatch(
        structure(list(mean = mu - lo, var = var), class = "MeanVarSpec"))
      return(shift_dist(inner, lo))
    }
    if (requested == "bounded") {
      # Truncation of a positive family — generic.
      return(solve_truncated_positive(name, mu, var, lo, hi))
    }
    if (requested == "half_left") {
      inner <- .dist_handlers[[name]]$dispatch(
        structure(list(mean = hi - mu, var = var), class = "MeanVarSpec"))
      return(flip_dist(inner, hi))
    }
  }

  if (natural == "unit") {
    if (requested == "bounded") {
      w <- hi - lo
      mu_std <- (mu - lo) / w
      var_std <- var / w^2
      inner <- .dist_handlers[[name]]$dispatch(
        structure(list(mean = mu_std, var = var_std), class = "MeanVarSpec"))
      return(scale_dist(inner, lo, w))
    }
  }

  stop(sprintf("Unsupported support combination: %s on (%g, %g)", name, lo, hi))
}


# Affine shift / flip / scale wrappers for an existing distsfactory_dist.

shift_dist <- function(inner, a) {
  out <- list(
    name = inner$name,
    params = c(inner$params, list(.shift = a)),
    d = function(x, log = FALSE) inner$d(x - a, log = log),
    p = function(q, lower.tail = TRUE, log.p = FALSE)
          inner$p(q - a, lower.tail = lower.tail, log.p = log.p),
    q = function(p, lower.tail = TRUE, log.p = FALSE)
          inner$q(p, lower.tail = lower.tail, log.p = log.p) + a,
    r = function(n) inner$r(n) + a
  )
  class(out) <- "distsfactory_dist"
  out
}

flip_dist <- function(inner, b) {
  out <- list(
    name = inner$name,
    params = c(inner$params, list(.flip = b)),
    # X = b - Y;  f_X(x) = f_Y(b - x);  F_X(x) = 1 - F_Y(b - x)
    d = function(x, log = FALSE) inner$d(b - x, log = log),
    p = function(q, lower.tail = TRUE, log.p = FALSE) {
      p <- 1 - inner$p(b - q)
      if (!lower.tail) p <- 1 - p
      if (log.p) log(p) else p
    },
    q = function(p, lower.tail = TRUE, log.p = FALSE) {
      if (log.p) p <- exp(p)
      if (!lower.tail) p <- 1 - p
      b - inner$q(1 - p)
    },
    r = function(n) b - inner$r(n)
  )
  class(out) <- "distsfactory_dist"
  out
}

scale_dist <- function(inner, a, w) {
  out <- list(
    name = inner$name,
    params = c(inner$params, list(.loc = a, .scale = w)),
    # X = a + w*Y;  f_X(x) = f_Y((x-a)/w) / w;  F_X(x) = F_Y((x-a)/w)
    d = function(x, log = FALSE) {
      if (log) inner$d((x - a) / w, log = TRUE) - log(w)
      else inner$d((x - a) / w) / w
    },
    p = function(q, lower.tail = TRUE, log.p = FALSE)
          inner$p((q - a) / w, lower.tail = lower.tail, log.p = log.p),
    q = function(p, lower.tail = TRUE, log.p = FALSE)
          a + w * inner$q(p, lower.tail = lower.tail, log.p = log.p),
    r = function(n) a + w * inner$r(n)
  )
  class(out) <- "distsfactory_dist"
  out
}


# ---------------------------------------------------------------------------
# Numerical solvers for truncated families.
# ---------------------------------------------------------------------------

# Underlying-family d/p/q/r function lookup for each locscale family. Used by
# the truncated-locscale solver to (re)build the underlying frozen distribution
# at trial (mu, sigma) values.
.locscale_dpq <- function(name) {
  switch(name,
    normal   = list(d = dnorm,   p = pnorm,   q = qnorm,   r = rnorm,
                    loc = "mean", scale = "sd",
                    var_of = function(p) p$sd^2,
                    mean_of = function(p) p$mean,
                    set = function(loc, sc) list(mean = loc, sd = sc)),
    laplace  = list(d = dlaplace_, p = plaplace_, q = qlaplace_, r = rlaplace_,
                    var_of = function(p) 2 * p$scale^2,
                    mean_of = function(p) p$location,
                    set = function(loc, sc) list(location = loc, scale = sc)),
    logistic = list(d = dlogis,  p = plogis,  q = qlogis,  r = rlogis,
                    var_of = function(p) (pi^2 / 3) * p$scale^2,
                    mean_of = function(p) p$location,
                    set = function(loc, sc) list(location = loc, scale = sc)),
    stop("Unknown locscale family: ", name)
  )
}

# Compute the (mean, var) of the truncated-locscale distribution with parent
# (loc, scale) on [lo, hi]. Done by numerical integration of x*f(x) and x^2*f(x)
# divided by the truncation mass.
.truncated_moments <- function(name, loc, scale, lo, hi) {
  fns <- .locscale_dpq(name)
  params <- fns$set(loc, scale)
  pdf <- function(x) do.call(fns$d, c(list(x), params))
  cdf <- function(x) do.call(fns$p, c(list(x), params))
  Z <- (if (is.infinite(hi) && hi > 0) 1 else cdf(hi)) -
       (if (is.infinite(lo) && lo < 0) 0 else cdf(lo))
  if (Z <= 0) return(c(NA, NA))
  # Integration: use a generous range capture. integrate() handles -Inf/Inf.
  ax <- if (is.infinite(lo)) -Inf else lo
  bx <- if (is.infinite(hi)) Inf else hi
  m  <- integrate(function(x) x * pdf(x) / Z, lower = ax, upper = bx,
                  rel.tol = 1e-10, abs.tol = 1e-12)$value
  m2 <- integrate(function(x) x^2 * pdf(x) / Z, lower = ax, upper = bx,
                  rel.tol = 1e-10, abs.tol = 1e-12)$value
  c(mean = m, var = m2 - m^2)
}

solve_truncated_locscale <- function(name, target_mean, target_var, lo, hi) {
  fns <- .locscale_dpq(name)
  # 2-D Newton on (loc, log_scale). Starting point: loc = midpoint of (lo, hi)
  # clamped, log_scale chosen from the half-width.
  loc0 <- if (is.finite(lo) && is.finite(hi)) (lo + hi) / 2
          else if (is.finite(lo)) target_mean
          else target_mean
  scale0 <- max(sqrt(target_var), 1e-3)
  start <- c(loc0, log(scale0))
  obj <- function(params) {
    loc <- params[1]; sc <- exp(params[2])
    mv <- .truncated_moments(name, loc, sc, lo, hi)
    c(mv["mean"] - target_mean, mv["var"] - target_var)
  }
  sol <- newton_2d(obj, start, maxiter = 200, tol = 1e-9, h = 1e-5)
  loc <- sol[1]; sc <- exp(sol[2])
  inner_params <- fns$set(loc, sc)
  make_truncated(name, inner_params, fns$d, fns$p, fns$q, fns$r, lo, hi)
}

# Generic positive-family truncation: assume base R d/p/q/r for the family,
# look up via .truncated_positive_funcs.
.positive_funcs <- function(name) {
  switch(name,
    gamma       = list(d = dgamma,  p = pgamma,  q = qgamma,  r = rgamma,
                       from_mv = gamma_from_mean_var,
                       params_of = function(d) list(shape = d$params$shape, rate = d$params$rate)),
    exponential = list(d = dexp,    p = pexp,    q = qexp,    r = rexp,
                       from_mv = exponential_from_mean_var %||% function(...) NULL,
                       params_of = function(d) list(rate = d$params$rate)),
    lognormal   = list(d = dlnorm,  p = plnorm,  q = qlnorm,  r = rlnorm,
                       from_mv = lognormal_from_mean_var,
                       params_of = function(d) list(meanlog = d$params$meanlog, sdlog = d$params$sdlog)),
    weibull     = list(d = dweibull, p = pweibull, q = qweibull, r = rweibull,
                       from_mv = weibull_from_mean_var,
                       params_of = function(d) list(shape = d$params$shape, scale = d$params$scale)),
    stop("Truncation of ", name, " on bounded support not implemented")
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

solve_truncated_positive <- function(name, target_mean, target_var, lo, hi) {
  fns <- .positive_funcs(name)
  # 2-D Newton on (log_param1, log_param2) — for two-parameter families.
  # Use the from_mean_var solution as a starting point (untruncated).
  start_dist <- fns$from_mv(target_mean, target_var)
  inner_params <- fns$params_of(start_dist)
  # Simple iteration: refine by adjusting the underlying mean and var so the
  # *truncated* moments match the target. Use Newton in (loc, log_scale)-like
  # space — but each family has its own parametrization, so we solve in the
  # natural (mean, var) of the *untruncated* family.
  obj <- function(params) {
    m_un <- params[1]; v_un <- exp(params[2])
    inner_dist <- fns$from_mv(m_un, v_un)
    ip <- fns$params_of(inner_dist)
    pdf <- function(x) do.call(fns$d, c(list(x), ip))
    cdf <- function(x) do.call(fns$p, c(list(x), ip))
    Z <- (if (is.infinite(hi)) 1 else cdf(hi)) -
         (if (lo <= 0) 0 else cdf(lo))
    ax <- max(lo, 0)
    bx <- if (is.infinite(hi)) Inf else hi
    m  <- integrate(function(x) x * pdf(x) / Z, lower = ax, upper = bx,
                    rel.tol = 1e-9)$value
    m2 <- integrate(function(x) x^2 * pdf(x) / Z, lower = ax, upper = bx,
                    rel.tol = 1e-9)$value
    c(m - target_mean, (m2 - m^2) - target_var)
  }
  start <- c(target_mean, log(target_var))
  sol <- newton_2d(obj, start, maxiter = 200, tol = 1e-9, h = 1e-5)
  m_un <- sol[1]; v_un <- exp(sol[2])
  inner_dist <- fns$from_mv(m_un, v_un)
  ip <- fns$params_of(inner_dist)
  make_truncated(name, ip, fns$d, fns$p, fns$q, fns$r, lo, hi)
}
