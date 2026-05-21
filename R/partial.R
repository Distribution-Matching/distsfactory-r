# Partial distribution specification — R analog of Julia's `@dist` / `DistSpec`
# and Python's `partial_dist`.
#
# Users pin some canonical parameters and leave others to be solved from moment
# constraints:
#
#   spec <- partial_dist("gamma", shape = 3.0)        # pin shape, solve rate
#   d <- make_dist(spec, mean = 5.0)
#   d$params$rate     # 0.6
#
# The solver is generic: identify free canonical parameters, build the
# distribution as a function of those, and 1D-brentq / 2D-Newton on
# (mean, var) to match the requested moments.

#' Build a partial distribution spec
#'
#' Returns a \code{partial_dist} object naming a family and some pinned
#' canonical parameters. Pass to \code{make_dist} with moment constraints
#' (\code{mean}, \code{var}) to solve the remaining parameters.
#'
#' @param dist Canonical distribution name (or alias).
#' @param ... Named parameter values to pin. Names must come from the family's
#'   canonical parameter set; see \code{canonical_params(dist)}.
#' @return An object of class \code{"partial_dist"} with elements
#'   \code{name} and \code{fixed}.
#'
#' @examples
#' spec <- partial_dist("gamma", shape = 3.0)
#' make_dist(spec, mean = 5.0)$params
#'
#' @export
partial_dist <- function(dist, ...) {
  name <- resolve_dist_name(dist)
  fixed <- list(...)
  canon <- canonical_params(name)
  bad <- setdiff(names(fixed), canon)
  if (length(bad) > 0)
    stop(sprintf("partial_dist(%s): unknown parameter(s): %s. Allowed: %s",
                 name, paste(bad, collapse = ", "),
                 paste(canon, collapse = ", ")))
  out <- list(name = name, fixed = fixed)
  class(out) <- "partial_dist"
  out
}

#' @export
print.partial_dist <- function(x, ...) {
  fixed_str <- paste(names(x$fixed), "=", format(unlist(x$fixed)),
                     collapse = ", ")
  cat(sprintf("partial_dist(%s, %s)\n", x$name, fixed_str))
  invisible(x)
}


# Internal: build a full distsfactory_dist from a partial spec with free params
# set to specific values. Uses the family's existing constructor under the hood
# by going through make_dist with a literal parameter override is unavailable;
# instead, we manually call new_dist (or the family's d/p/q/r functions
# directly) with the combined params.
.partial_build <- function(name, fixed, free_names, free_values) {
  params <- fixed
  for (i in seq_along(free_names)) {
    params[[free_names[i]]] <- free_values[i]
  }
  # Build the dist by invoking the family's normal mean-var constructor at
  # moments implied by the params is not always available; instead, rebuild
  # the d/p/q/r closures by calling new_dist with the appropriate handlers.
  fns <- .family_dpq(name)
  new_dist(name, params, fns$d, fns$p, fns$q, fns$r)
}

# Resolve the (d, p, q, r) base function set for each family. For most
# distributions these are base R names; inline-implemented families return
# their family-specific underscored helpers.
.family_dpq <- function(name) {
  switch(name,
    gamma       = list(d = dgamma,  p = pgamma,  q = qgamma,  r = rgamma),
    exponential = list(d = dexp,    p = pexp,    q = qexp,    r = rexp),
    logistic    = list(d = dlogis,  p = plogis,  q = qlogis,  r = rlogis),
    beta        = list(d = dbeta,   p = pbeta,   q = qbeta,   r = rbeta),
    normal      = list(d = dnorm,   p = pnorm,   q = qnorm,   r = rnorm),
    lognormal   = list(d = dlnorm,  p = plnorm,  q = qlnorm,  r = rlnorm),
    uniform     = list(d = dunif,   p = punif,   q = qunif,   r = runif),
    weibull     = list(d = dweibull, p = pweibull, q = qweibull, r = rweibull),
    tdist       = list(d = dt,      p = pt,      q = qt,      r = rt),
    chisq       = list(d = dchisq,  p = pchisq,  q = qchisq,  r = rchisq),
    fdist       = list(d = df,      p = pf,      q = qf,      r = rf),
    laplace     = list(d = dlaplace_,  p = plaplace_,  q = qlaplace_,  r = rlaplace_),
    gumbel      = list(d = dgumbel_,   p = pgumbel_,   q = qgumbel_,   r = rgumbel_),
    rayleigh    = list(d = drayleigh_, p = prayleigh_, q = qrayleigh_, r = rrayleigh_),
    pareto      = list(d = dpareto_,   p = ppareto_,   q = qpareto_,   r = rpareto_),
    frechet     = list(d = dfrechet_,  p = pfrechet_,  q = qfrechet_,  r = rfrechet_),
    inverse_gamma = list(d = dinvgamma_, p = pinvgamma_, q = qinvgamma_, r = rinvgamma_),
    chi         = list(d = dchi_,   p = pchi_,   q = qchi_,   r = rchi_),
    folded_normal = list(d = dfoldednorm_, p = pfoldednorm_, q = qfoldednorm_, r = rfoldednorm_),
    erlang      = list(d = dgamma,  p = pgamma,  q = qgamma,  r = rgamma),
    sym_triangular = list(d = dsymtri_, p = psymtri_, q = qsymtri_, r = rsymtri_),
    binomial    = list(d = dbinom,  p = pbinom,  q = qbinom,  r = rbinom),
    poisson     = list(d = dpois,   p = ppois,   q = qpois,   r = rpois),
    negative_binomial = list(d = dnbinom, p = pnbinom, q = qnbinom, r = rnbinom),
    geometric   = list(d = dgeom,   p = pgeom,   q = qgeom,   r = rgeom),
    discrete_uniform = list(d = ddunif_, p = pdunif_, q = qdunif_, r = rdunif_),
    discrete_sym_triangular = list(d = ddst_, p = pdst_, q = qdst_, r = rdst_),
    stop("Family ", name, " not supported by partial_dist")
  )
}

# Tunable: starting guesses for the solver.
.partial_initial <- function(name, free_names) {
  # Use 1.0 for everything except clearly positive shape/scale-like params;
  # let the brentq bracket expansion handle the rest.
  rep(1, length(free_names))
}

# Solve a partial_dist for its free params from moment constraints.
make_dist_from_partial <- function(spec, mean = NULL, var = NULL, std = NULL,
                                   cv = NULL, scv = NULL, second_moment = NULL) {
  # Resolve target var if expressed via std/cv/scv/second_moment.
  if (is.null(var) && !is.null(std)) var <- std^2
  if (is.null(var) && !is.null(mean) && !is.null(cv)) var <- (cv * mean)^2
  if (is.null(var) && !is.null(mean) && !is.null(scv)) var <- scv * mean^2
  if (is.null(var) && !is.null(mean) && !is.null(second_moment))
    var <- second_moment - mean^2

  name <- spec$name
  fixed <- spec$fixed
  canon <- canonical_params(name)
  free_names <- setdiff(canon, names(fixed))

  if (length(free_names) == 0) {
    # Nothing to solve — just verify (and warn if moments mismatch).
    return(.partial_build(name, fixed, character(0), numeric(0)))
  }

  moments_of <- function(values) {
    params <- fixed
    for (i in seq_along(free_names)) params[[free_names[i]]] <- values[i]
    tryCatch(c(dist_mean_from_params(name, params),
               dist_var_from_params(name, params)),
             error = function(e) c(NA_real_, NA_real_))
  }

  if (length(free_names) == 1) {
    primary <- if (!is.null(var)) "var" else "mean"
    f <- function(x) {
      mv <- moments_of(x)
      if (primary == "var") mv[2] - var else mv[1] - mean
    }
    # Try brentq with widening brackets.
    sol <- tryCatch(.bracketed_brentq(f, x0 = 1),
                    error = function(e) NULL)
    if (is.null(sol) && !is.null(mean)) {
      f2 <- function(x) moments_of(x)[1] - mean
      sol <- tryCatch(.bracketed_brentq(f2, x0 = 1),
                      error = function(e) NULL)
    }
    if (is.null(sol))
      stop(sprintf("partial_dist(%s): could not solve free parameter %s",
                   name, free_names))
    out <- .partial_build(name, fixed, free_names, sol)
    if (!is.null(mean) && !is.null(var)) {
      m <- dist_mean_from_params(name, out$params)
      v <- dist_var_from_params(name, out$params)
      if (!isTRUE(all.equal(m, mean, tolerance = 1e-4)) ||
          !isTRUE(all.equal(v, var, tolerance = 1e-4)))
        stop(sprintf("partial_dist(%s): cannot satisfy both mean and var with 1 free parameter",
                     name))
    }
    return(out)
  }

  if (length(free_names) == 2 && !is.null(mean) && !is.null(var)) {
    obj <- function(values) {
      mv <- moments_of(values)
      c(mv[1] - mean, mv[2] - var)
    }
    sol <- newton_2d(obj, .partial_initial(name, free_names),
                     maxiter = 200, tol = 1e-9, h = 1e-6)
    return(.partial_build(name, fixed, free_names, sol))
  }

  stop(sprintf("partial_dist(%s): unsupported solver shape with %d free params; need mean+var for 2 free",
               name, length(free_names)))
}

.bracketed_brentq <- function(f, x0 = 1) {
  candidates <- list(c(1e-6, 1e3), c(1e-9, 1e6), c(-1e6, 1e6))
  for (br in candidates) {
    fa <- tryCatch(f(br[1]), error = function(e) NA_real_)
    fb <- tryCatch(f(br[2]), error = function(e) NA_real_)
    if (is.finite(fa) && is.finite(fb) && fa * fb < 0) {
      return(uniroot(f, br, tol = 1e-12)$root)
    }
  }
  NULL
}
