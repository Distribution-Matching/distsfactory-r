# Regression-protect every code example in the README.
# Each block is run verbatim (or close to it); only the assertions are added
# so a future README rewrite that breaks an example surfaces immediately.

test_that("README Quick start: gamma mean+var basic", {
  d <- make_dist("gamma", mean = 5, var = 3)
  expect_s3_class(d, "distsfactory_dist")
  expect_true(is.numeric(d$d(2)))
  expect_true(d$p(0.95) > 0 && d$p(0.95) < 1)
  expect_true(is.numeric(d$q(0.5)))
  expect_length(d$r(100), 100)
  expect_named(d$params, c("shape", "rate"))
})

test_that("README Quick start: truncated Normal on [-1, 4]", {
  d <- make_dist("normal", mean = 1.0, std = 0.8, support = c(-1, 4))
  expect_equal(d$name, "truncated_normal")
  expect_equal(d$mean(), 1.0, tolerance = 1e-6)
  expect_equal(d$std(),  0.8, tolerance = 1e-6)
  expect_equal(d$support, c(-1, 4))
  expect_equal(d$parent$params$mean, 0.9822, tolerance = 1e-3)
  expect_equal(d$parent$params$sd,   0.8232, tolerance = 1e-3)
  expect_equal(d$p(-1), 0, tolerance = 1e-9)
  expect_equal(d$p(4),  1, tolerance = 1e-9)
})

# --- Specification styles block ----------------------------------------

test_that("README spec styles: moment alternatives all run", {
  expect_silent(make_dist("gamma", mean = 5, var = 3))
  expect_silent(make_dist("gamma", mean = 5, std = 2))
  expect_silent(make_dist("gamma", mean = 5, cv = 0.5))
  expect_silent(make_dist("gamma", mean = 4, scv = 0.5))
  expect_silent(make_dist("gamma", mean = 5, second_moment = 28))
  expect_silent(make_dist("exponential", mean = 3))
})

test_that("README spec styles: quantile-based all run", {
  expect_silent(make_dist("exponential", median = 2.0))
  expect_silent(make_dist("logistic", q1 = 2, q3 = 8))
  expect_silent(make_dist("normal", q1 = -1, q3 = 1))
  expect_silent(make_dist("gamma", quantiles = list(c(0.1, 1.0), c(0.9, 10.0))))
  expect_silent(make_dist("beta", mean = 0.4, median = 0.38))
  expect_silent(make_dist("normal", median = 5, iqr = 2))
})

test_that("README spec styles: mode-based all run", {
  expect_silent(make_dist("rayleigh", mode = 2))
  expect_silent(make_dist("gamma", mean = 5, mode = 3))
  expect_silent(make_dist("beta", mean = 0.4, mode = 0.35))
  expect_silent(make_dist("gamma", mode = 3, iqr = 4))
  expect_silent(make_dist("normal", mode = 3, var = 4))
  expect_silent(make_dist("logistic", mode = 5, iqr = 4))
})

test_that("README spec styles: 3-param triangular", {
  expect_silent(make_dist("triangular", mean = 5, var = 2, mode = 4))
  expect_silent(make_dist("discrete_triangular", mean = 5, var = 2, mode = 5))
})

# --- Support block -----------------------------------------------------

test_that("README support: affine shift gamma [3, Inf)", {
  d <- make_dist("gamma", mean = 8, var = 3, support = c(3, Inf))
  expect_equal(d$p(3), 0, tolerance = 1e-9)
})

test_that("README support: affine flip gamma (-Inf, 10]", {
  d <- make_dist("gamma", mean = 5, var = 3, support = c(-Inf, 10))
  expect_equal(d$p(10), 1, tolerance = 1e-9)
})

test_that("README support: affine scale beta on [2, 7]", {
  d <- make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))
  expect_equal(d$p(2), 0, tolerance = 1e-9)
  expect_equal(d$p(7), 1, tolerance = 1e-9)
})

test_that("README support: truncation Normal on [-0.5, 0.5]", {
  d <- make_dist("normal", mean = 0.1, var = 0.05, support = c(-0.5, 0.5))
  expect_equal(d$p(-0.5), 0, tolerance = 1e-9)
  expect_equal(d$p(0.5),  1, tolerance = 1e-9)
})

test_that("README support: truncation gamma on [0, 10]", {
  expect_silent(make_dist("gamma", mean = 3, var = 1, support = c(0, 10)))
})

test_that("README support: truncation beta on [0.2, 0.8]", {
  expect_silent(make_dist("beta",  mean = 0.5, var = 0.02, support = c(0.2, 0.8)))
})

# --- partial_dist block ------------------------------------------------

test_that("README partial: pin gamma shape, solve rate from mean", {
  spec <- partial_dist("gamma", shape = 3.0)
  d <- make_dist(spec, mean = 5.0)
  expect_equal(d$params$shape, 3.0)
  expect_equal(d$params$rate, 0.6, tolerance = 1e-9)
})

test_that("README partial: pin gamma shape, solve from variance", {
  spec <- partial_dist("gamma", shape = 3.0)
  d <- make_dist(spec, var = 3.0)
  expect_equal(d$params$shape, 3.0)
  expect_equal(d$params$shape / d$params$rate^2, 3, tolerance = 1e-9)
})

test_that("README partial: pin normal mean, solve sd from var", {
  spec <- partial_dist("normal", mean = 2.0)
  d <- make_dist(spec, var = 4.0)
  expect_equal(d$params$mean, 2.0)
  expect_equal(d$params$sd, 2.0, tolerance = 1e-9)
})

test_that("README partial: no pins — solve both from mean+var", {
  d <- make_dist(partial_dist("gamma"), mean = 5, var = 3)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-6)
  expect_equal(d$params$shape / d$params$rate^2, 3, tolerance = 1e-6)
})

# --- Feasibility block -------------------------------------------------

test_that("README feasibility: beta", {
  expect_true(dist_exists("beta", mean = 0.5, var = 0.1))
  expect_false(dist_exists("beta", mean = 0.5, var = 0.3))
})

test_that("README feasibility: exponential", {
  expect_true(dist_exists("exponential", mean = 2.5, var = 6.25))
  expect_false(dist_exists("exponential", mean = 2.5, var = 1.5))
})

test_that("README feasibility: tdist requires mean=0", {
  expect_false(dist_exists("tdist", mean = 1, var = 2))
})

# --- Discovery block ---------------------------------------------------

test_that("README discovery: mean+var enumeration", {
  dists <- available_distributions(mean = 5, var = 3)
  expect_true("logistic" %in% dists)
  expect_true("gamma" %in% dists)
})

test_that("README discovery: support-filtered enumeration", {
  dists <- available_distributions(mean = 3.5, var = 0.5, support = c(2, 7))
  expect_true("beta" %in% dists)
})
