# Langevin feasibility envelope for truncated location-scale families.
#
# Port of src/langevin.jl in DistributionsFactories.jl. For any continuous,
# unimodal, location-scale family with exponential tails (Normal, Laplace,
# Logistic) truncated to a bounded interval [a, b], the set of achievable
# (mean, variance) pairs forms a dome whose upper envelope coincides with the
# truncated-exponential family. The maximum variance at mean mu on [a, b] is
#
#     sigma2_max(mu) = w^2 * L'(L^{-1}((c - mu) / w))
#
# where c = (a + b) / 2, w = (b - a) / 2, and L is the Langevin function
# L(x) = coth(x) - 1/x.
#
# For half-truncated cases (lo finite, hi = +Inf), the exponential-tail bound
# gives var < (mu - lo)^2; the Laplace boundary is attained when the parent
# location is past the truncation point.

langevin <- function(x) {
  if (abs(x) < 1e-4)
    return(x / 3 - x^3 / 45 + 2 * x^5 / 945 - x^7 / 4725)
  1 / tanh(x) - 1 / x
}

langevin_deriv <- function(x) {
  if (abs(x) < 1e-3)
    return(1 / 3 - x^2 / 15 + 2 * x^4 / 189 - x^6 / 675)
  s <- sinh(x)
  1 / x^2 - 1 / (s * s)
}

inv_langevin <- function(y) {
  if (!(abs(y) < 1))
    stop(sprintf("inv_langevin: |y| must be < 1 (got %g)", y))
  if (y == 0) return(0)
  # Cohen's (1991) [3/2] Pade initial guess.
  z <- y * (3 - y^2) / (1 - y^2)
  for (i in seq_len(50)) {
    f <- langevin(z) - y
    if (abs(f) < 1e-14) return(z)
    z <- z - f / langevin_deriv(z)
  }
  z
}

truncexp_max_var <- function(a, b, mu) {
  if (!(a < b)) stop(sprintf("truncexp_max_var: require a < b (got a=%g, b=%g)", a, b))
  if (!(a < mu && mu < b))
    stop(sprintf("truncexp_max_var: mu must lie strictly in (%g, %g) (got %g)", a, b, mu))
  c <- (a + b) / 2
  w <- (b - a) / 2
  z <- inv_langevin((c - mu) / w)
  w^2 * langevin_deriv(z)
}

# Feasibility predicate for truncated Normal / Laplace / Logistic on [lo, hi].
# Returns TRUE iff a distribution in the named family with the given moments
# can be placed on [lo, hi] (allowing lo = -Inf or hi = +Inf for half-truncated
# cases). The full bound is the Langevin dome on bounded intervals; on
# half-infinite intervals it reduces to var < (mu - lo)^2 (or var < (hi - mu)^2).
truncated_locscale_exists <- function(name, mu, var, lo, hi) {
  if (!is.finite(mu) || !is.finite(var) || var <= 0) return(FALSE)
  if (!(lo < hi)) return(FALSE)
  if (!(lo < mu && mu < hi)) return(FALSE)

  if (is.finite(lo) && is.finite(hi)) {
    bound <- truncexp_max_var(lo, hi, mu)
    return(var < bound)
  }
  if (is.finite(lo) && is.infinite(hi)) {
    bound <- (mu - lo)^2
    if (name == "laplace") return(var <= bound)   # boundary attained
    return(var < bound)
  }
  if (is.infinite(lo) && is.finite(hi)) {
    bound <- (hi - mu)^2
    if (name == "laplace") return(var <= bound)
    return(var < bound)
  }
  # both infinite -> reduces to untruncated family
  TRUE
}
