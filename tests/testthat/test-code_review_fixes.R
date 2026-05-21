# Regression tests for the issues caught in the package-wide code review.

# --- 1. support= wrappers honour log / lower.tail / log.p -----------------

test_that("truncated d/p/q honor log, lower.tail, log.p", {
  d <- make_dist("normal", mean = 0, var = 1, support = c(-2, 2))
  # density
  expect_equal(d$d(0, log = TRUE), log(d$d(0)), tolerance = 1e-12)
  # cdf
  expect_equal(d$p(0, lower.tail = FALSE), 1 - d$p(0), tolerance = 1e-12)
  expect_equal(d$p(0, log.p = TRUE), log(d$p(0)), tolerance = 1e-12)
  # qf
  expect_equal(d$q(log(0.3), log.p = TRUE), d$q(0.3), tolerance = 1e-9)
  expect_equal(d$q(0.7, lower.tail = FALSE), d$q(0.3), tolerance = 1e-9)
})

test_that("shifted dist (positive on half-right) honors log/lower.tail/log.p", {
  d <- make_dist("gamma", mean = 8, var = 4, support = c(3, Inf))
  x_in <- 5
  expect_equal(d$d(x_in, log = TRUE), log(d$d(x_in)), tolerance = 1e-12)
  expect_equal(d$p(x_in, lower.tail = FALSE), 1 - d$p(x_in), tolerance = 1e-12)
  expect_equal(d$p(x_in, log.p = TRUE), log(d$p(x_in)), tolerance = 1e-12)
  expect_equal(d$q(0.3, lower.tail = FALSE), d$q(0.7), tolerance = 1e-9)
})

test_that("flipped dist on (-Inf, hi] honors flags", {
  d <- make_dist("gamma", mean = 5, var = 3, support = c(-Inf, 10))
  expect_equal(d$d(8, log = TRUE), log(d$d(8)), tolerance = 1e-12)
  expect_equal(d$p(8, lower.tail = FALSE), 1 - d$p(8), tolerance = 1e-12)
  expect_equal(d$p(8, log.p = TRUE), log(d$p(8)), tolerance = 1e-12)
})

test_that("scaled dist (beta on bounded) honors flags", {
  d <- make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))
  # 1/w factor in log space
  expect_equal(d$d(3.5, log = TRUE), log(d$d(3.5)), tolerance = 1e-12)
  expect_equal(d$p(3.5, lower.tail = FALSE), 1 - d$p(3.5), tolerance = 1e-12)
})


# --- 2. Conflicting dispersion measures are caught ------------------------

test_that("conflicting var and std are rejected", {
  expect_error(make_dist("gamma", mean = 5, var = 3, std = 2),
               "Conflicting dispersion measures")
})

test_that("consistent var and std are accepted (no error)", {
  expect_silent(make_dist("gamma", mean = 5, var = 4, std = 2))
})

test_that("conflicting cv and var rejected", {
  # var = (cv*mean)^2 = (0.5*5)^2 = 6.25; passing var = 3 conflicts.
  expect_error(make_dist("gamma", mean = 5, var = 3, cv = 0.5),
               "Conflicting dispersion measures")
})

test_that("conflicting second_moment and var rejected", {
  # second_moment - mean^2 = 28 - 25 = 3, so var = 4 conflicts.
  expect_error(make_dist("gamma", mean = 5, var = 4, second_moment = 28),
               "Conflicting dispersion measures")
})


# --- 3. Non-finite moments are rejected -----------------------------------

test_that("NA mean is rejected", {
  expect_error(make_dist("normal", mean = NA, var = 1), "must be finite")
})

test_that("NA var is rejected", {
  expect_error(make_dist("normal", mean = 0, var = NA), "must be finite")
})

test_that("NaN std is rejected", {
  expect_error(make_dist("normal", mean = 0, std = NaN), "must be finite")
})

test_that("Inf quantile is rejected", {
  expect_error(make_dist("normal", q1 = -Inf, q3 = 1), "must be finite")
})


# --- 4. print method handles vector params --------------------------------

test_that("print handles vector parameters without flattening", {
  # Craft a custom dist with a vector param to verify the format.
  obj <- list(
    name = "fake",
    params = list(a = c(1, 2, 3), b = 5),
    d = function(x) x, p = function(q) q, q = function(p) p, r = function(n) n
  )
  class(obj) <- "distsfactory_dist"
  out <- capture.output(print(obj))
  expect_match(out, "a = c\\(1, 2, 3\\)")
  expect_match(out, "b = 5")
})

test_that("print still works for scalar params (regression)", {
  d <- make_dist("gamma", mean = 5, var = 3)
  out <- capture.output(print(d))
  expect_match(out, "distsfactory: gamma")
  expect_match(out, "shape")
  expect_match(out, "rate")
})


# --- 5. available_distributions with no args lists registry --------------

test_that("available_distributions() with no args returns the registry", {
  dists <- available_distributions()
  expect_true(length(dists) >= 25)
  expect_true("gamma" %in% dists)
  expect_true("normal" %in% dists)
  expect_true("beta" %in% dists)
  expect_true("poisson" %in% dists)
})


# --- 6. partial_dist now knows about triangular + discrete_triangular ---

test_that("partial_dist on triangular with pinned mode + 2 free (a, b)", {
  # mean = 5, var = 2, mode pinned via partial_dist. 2 free params (a, b).
  spec <- partial_dist("triangular", c = 4)
  d <- make_dist(spec, mean = 5, var = 2)
  expect_equal(d$params$c, 4)
  expect_equal((d$params$a + d$params$b + d$params$c) / 3, 5, tolerance = 1e-5)
  achieved_var <- (d$params$a^2 + d$params$b^2 + d$params$c^2 -
                   d$params$a * d$params$b - d$params$a * d$params$c -
                   d$params$b * d$params$c) / 18
  expect_equal(achieved_var, 2, tolerance = 1e-5)
})


# --- 7. Pareto mean+var still works after cleanup ------------------------

test_that("Pareto mean+var (post-cleanup) computes correct alpha", {
  d <- make_dist("pareto", mean = 5, var = 3)
  a <- d$params$shape; xm <- d$params$scale
  # Verify formulas: mean = a*xm/(a-1), var = xm^2*a / ((a-1)^2*(a-2))
  expect_equal(a * xm / (a - 1), 5, tolerance = 1e-9)
  expect_equal(xm^2 * a / ((a - 1)^2 * (a - 2)), 3, tolerance = 1e-9)
  expect_gt(a, 2)
})


# --- 8. Bracket improvements: partial_dist on pareto -----------------

test_that("partial_dist(pareto, pin shape > 2) solves cleanly", {
  # mean = a*xm/(a-1). Pin a=3, target mean=6 -> xm = 4.
  spec <- partial_dist("pareto", shape = 3)
  d <- make_dist(spec, mean = 6)
  expect_equal(d$params$shape, 3)
  expect_equal(d$params$scale, 4, tolerance = 1e-6)
})

test_that("partial_dist(frechet, pin shape > 2) solves cleanly", {
  # mean = scale * Gamma(1 - 1/alpha). alpha=3 -> mean = scale * Gamma(2/3).
  spec <- partial_dist("frechet", shape = 3)
  d <- make_dist(spec, mean = 5)
  expect_equal(d$params$shape, 3)
  expect_equal(d$params$scale * gamma(1 - 1/3), 5, tolerance = 1e-6)
})

test_that("partial_dist(inverse_gamma, pin shape > 2) solves cleanly", {
  spec <- partial_dist("inverse_gamma", shape = 3)
  d <- make_dist(spec, mean = 5)
  expect_equal(d$params$shape, 3)
  expect_equal(d$params$scale / (d$params$shape - 1), 5, tolerance = 1e-9)
})
