test_that("partial_dist pinning gamma shape solves rate from mean", {
  spec <- partial_dist("gamma", shape = 3.0)
  d <- make_dist(spec, mean = 5.0)
  # mean = shape / rate -> rate = shape / mean = 3/5 = 0.6
  expect_equal(d$params$shape, 3.0)
  expect_equal(d$params$rate, 0.6, tolerance = 1e-9)
})

test_that("partial_dist pinning normal mean solves sd from var", {
  spec <- partial_dist("normal", mean = 1.0)
  d <- make_dist(spec, var = 4.0)
  expect_equal(d$params$mean, 1.0)
  expect_equal(d$params$sd, 2.0, tolerance = 1e-9)
})

test_that("partial_dist with no fixed params (2 free) solves mean+var", {
  spec <- partial_dist("gamma")  # nothing pinned
  d <- make_dist(spec, mean = 5, var = 3)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-6)
  expect_equal(d$params$shape / d$params$rate^2, 3, tolerance = 1e-6)
})

test_that("partial_dist normal alias / unknown param raises", {
  expect_error(partial_dist("gamma", nonsense = 1), "unknown parameter")
})

test_that("partial_dist with 1-param family + var-only solves", {
  spec <- partial_dist("exponential")  # rate is free, single param
  d <- make_dist(spec, mean = 5)
  expect_equal(d$params$rate, 0.2, tolerance = 1e-9)
})

test_that("partial_dist returns distsfactory_dist", {
  spec <- partial_dist("gamma", shape = 2)
  d <- make_dist(spec, mean = 4)
  expect_s3_class(d, "distsfactory_dist")
  expect_true(is.function(d$d))
  expect_true(is.function(d$p))
})

test_that("partial_dist print method", {
  spec <- partial_dist("gamma", shape = 2)
  expect_output(print(spec), "partial_dist\\(gamma")
})


# ---------------------------------------------------------------------------
# Coverage sweep: one happy-path test per family with at least one canonical
# parameter. Each case builds a partial_dist (some pinned, some free), solves
# from moment constraints, and verifies that the resolved parameters produce
# moments that match the requested targets to a generous tolerance. For
# 2-parameter families both the "pin one, solve the other" and "pin none,
# solve both" patterns are exercised where they converge with the default
# starting guess (rep(1, k)).
# ---------------------------------------------------------------------------

# Helper: build via partial_dist, then verify achieved (mean, var) against
# the requested targets.
.check_partial <- function(name, fixed = list(), mean = NULL, var = NULL,
                           tol = 1e-5) {
  spec <- do.call(partial_dist, c(list(name), fixed))
  d <- make_dist(spec, mean = mean, var = var)
  achieved_mean <- dist_mean_from_params(name, d$params)
  achieved_var  <- dist_var_from_params(name, d$params)
  if (!is.null(mean))
    testthat::expect_equal(achieved_mean, mean, tolerance = tol,
                           info = sprintf("%s: achieved mean", name))
  if (!is.null(var))
    testthat::expect_equal(achieved_var, var, tolerance = tol,
                           info = sprintf("%s: achieved var", name))
  d
}

# --- 1-parameter continuous ----------------------------------------------

test_that("partial_dist coverage: exponential (1 free)", {
  .check_partial("exponential", mean = 4)
})

test_that("partial_dist coverage: tdist (1 free, var only)", {
  # TDist has mean = 0; var = df/(df-2) for df > 2. Pin nothing, solve from var.
  .check_partial("tdist", var = 3)
})

test_that("partial_dist coverage: chisq (1 free)", {
  .check_partial("chisq", mean = 6)
})

test_that("partial_dist coverage: chi (1 free)", {
  # E[Chi(df)] = sqrt(2)*Gamma((df+1)/2)/Gamma(df/2). Solve from mean.
  target <- sqrt(2) * exp(lgamma(2) - lgamma(1.5))   # E[Chi(3)] ≈ 1.5958
  .check_partial("chi", mean = target)
})

test_that("partial_dist coverage: rayleigh (1 free)", {
  .check_partial("rayleigh", mean = 5)
})

test_that("partial_dist coverage: poisson (1 free)", {
  .check_partial("poisson", mean = 7)
})

test_that("partial_dist coverage: geometric (1 free)", {
  .check_partial("geometric", mean = 3)
})

# --- 2-parameter continuous: pin one, solve other -----------------------

test_that("partial_dist coverage: gamma (pin shape)", {
  d <- .check_partial("gamma", fixed = list(shape = 3.0), mean = 5)
  expect_equal(d$params$shape, 3.0)
})

test_that("partial_dist coverage: gamma (pin rate)", {
  d <- .check_partial("gamma", fixed = list(rate = 0.5), mean = 6)
  expect_equal(d$params$rate, 0.5)
})

test_that("partial_dist coverage: normal (pin mean)", {
  d <- .check_partial("normal", fixed = list(mean = 2), var = 9)
  expect_equal(d$params$mean, 2)
  expect_equal(d$params$sd, 3, tolerance = 1e-9)
})

test_that("partial_dist coverage: normal (pin sd)", {
  d <- .check_partial("normal", fixed = list(sd = 2), mean = 4)
  expect_equal(d$params$sd, 2)
})

test_that("partial_dist coverage: logistic (pin location)", {
  .check_partial("logistic", fixed = list(location = 1), var = 10)
})

test_that("partial_dist coverage: laplace (pin location)", {
  .check_partial("laplace", fixed = list(location = 0), var = 2)
})

test_that("partial_dist coverage: beta (pin shape1)", {
  # mean = a/(a+b) -> b = a/mean - a = 2/0.4 - 2 = 3
  d <- .check_partial("beta", fixed = list(shape1 = 2), mean = 0.4)
  expect_equal(d$params$shape2, 3, tolerance = 1e-9)
})

test_that("partial_dist coverage: lognormal (pin meanlog)", {
  # meanlog = 0 -> mean = exp(sigma^2/2). Pick a sigma to back-solve.
  spec <- partial_dist("lognormal", meanlog = 0)
  d <- make_dist(spec, mean = 2)
  expect_equal(d$params$meanlog, 0)
  expect_equal(exp(d$params$meanlog + d$params$sdlog^2 / 2), 2,
               tolerance = 1e-9)
})

test_that("partial_dist coverage: uniform (pin min)", {
  d <- .check_partial("uniform", fixed = list(min = 0), mean = 5)
  expect_equal(d$params$min, 0)
  expect_equal(d$params$max, 10, tolerance = 1e-9)
})

test_that("partial_dist coverage: weibull (pin shape)", {
  # mean = scale * Gamma(1 + 1/k). With k=2, scale = 2*mean/sqrt(pi)
  d <- .check_partial("weibull", fixed = list(shape = 2), mean = 5)
  expect_equal(d$params$shape, 2)
  expect_equal(d$params$scale, 2 * 5 / sqrt(pi), tolerance = 1e-9)
})

test_that("partial_dist coverage: pareto (pin shape)", {
  # mean = shape*scale/(shape-1). shape=3 -> scale = 2*mean/3
  d <- .check_partial("pareto", fixed = list(shape = 3), mean = 6)
  expect_equal(d$params$shape, 3)
  expect_equal(d$params$scale, 2 * 6 / 3, tolerance = 1e-9)
})

test_that("partial_dist coverage: inverse_gamma (pin shape)", {
  # mean = scale/(shape-1). shape=3 -> scale = 2*mean
  d <- .check_partial("inverse_gamma", fixed = list(shape = 3), mean = 5)
  expect_equal(d$params$shape, 3)
  expect_equal(d$params$scale, 10, tolerance = 1e-9)
})

test_that("partial_dist coverage: gumbel (pin scale)", {
  # mean = location + scale * EulerGamma. scale=1 -> location = mean - gamma
  euler <- -digamma(1)
  d <- .check_partial("gumbel", fixed = list(scale = 1), mean = 5)
  expect_equal(d$params$location, 5 - euler, tolerance = 1e-9)
})

test_that("partial_dist coverage: sym_triangular (pin location)", {
  d <- .check_partial("sym_triangular", fixed = list(location = 0), var = 6)
  expect_equal(d$params$location, 0)
  expect_equal(d$params$scale, 6, tolerance = 1e-9)
})

# --- 2-parameter discrete: pin one, solve other -------------------------

test_that("partial_dist coverage: binomial (pin size)", {
  # mean = size*prob -> prob = mean/size = 5/10 = 0.5
  d <- .check_partial("binomial", fixed = list(size = 10), mean = 5)
  expect_equal(d$params$prob, 0.5, tolerance = 1e-9)
})

test_that("partial_dist coverage: negative_binomial (pin prob)", {
  # mean = size*(1-prob)/prob. prob=0.5 -> size = mean = 5
  d <- .check_partial("negative_binomial", fixed = list(prob = 0.5), mean = 5)
  expect_equal(d$params$prob, 0.5)
  expect_equal(d$params$size, 5, tolerance = 1e-9)
})

test_that("partial_dist coverage: discrete_uniform (pin min)", {
  # mean = (min+max)/2 -> max = 2*mean - min = 10
  d <- .check_partial("discrete_uniform", fixed = list(min = 0), mean = 5)
  expect_equal(d$params$min, 0)
  expect_equal(d$params$max, 10, tolerance = 1e-9)
})

test_that("partial_dist coverage: discrete_sym_triangular (pin mu)", {
  # var = n*(n+2)/6. var=4 -> n*(n+2)=24 -> n=4
  d <- .check_partial("discrete_sym_triangular", fixed = list(mu = 5), var = 4)
  expect_equal(d$params$mu, 5)
  expect_equal(d$params$n, 4, tolerance = 1e-6)
})

# --- 2-free Newton-on-(mean, var) ---------------------------------------

test_that("partial_dist coverage: 2-free Newton — gamma", {
  d <- .check_partial("gamma", mean = 5, var = 3, tol = 1e-5)
  expect_equal(d$params$shape, 25 / 3, tolerance = 1e-4)
  expect_equal(d$params$rate, 5 / 3,  tolerance = 1e-4)
})

test_that("partial_dist coverage: 2-free Newton — normal", {
  .check_partial("normal", mean = 0, var = 4)
})

test_that("partial_dist coverage: 2-free Newton — beta", {
  .check_partial("beta", mean = 0.5, var = 0.05)
})

test_that("partial_dist coverage: 2-free Newton — logistic", {
  .check_partial("logistic", mean = 3, var = 10)
})

test_that("partial_dist coverage: 2-free Newton — uniform", {
  .check_partial("uniform", mean = 5, var = 3)
})

test_that("partial_dist coverage: 2-free Newton — laplace", {
  .check_partial("laplace", mean = 5, var = 8)
})

# --- Unsupported / error paths ------------------------------------------

test_that("partial_dist with 3 free params errors cleanly", {
  spec <- partial_dist("triangular")
  expect_error(make_dist(spec, mean = 5, var = 2, mode = 4),
               "unsupported solver shape")
})

test_that("partial_dist with all params pinned skips solving", {
  spec <- partial_dist("gamma", shape = 2, rate = 0.5)
  d <- make_dist(spec)
  expect_equal(d$params$shape, 2)
  expect_equal(d$params$rate, 0.5)
})
