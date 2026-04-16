test_that("logistic from mean+var basic", {
  d <- make_dist("logistic", mean = 0, var = 3)
  expect_equal(d$params$location, 0, tolerance = 1e-6)
  # var = scale^2 * pi^2 / 3
  computed_var <- d$params$scale^2 * pi^2 / 3
  expect_equal(computed_var, 3, tolerance = 1e-6)
})

test_that("logistic from mean+var nonzero mean", {
  d <- make_dist("logistic", mean = 10, var = 5)
  expect_equal(d$params$location, 10, tolerance = 1e-6)
  computed_var <- d$params$scale^2 * pi^2 / 3
  expect_equal(computed_var, 5, tolerance = 1e-6)
})

test_that("logistic from mean+std", {
  d <- make_dist("logistic", mean = 0, std = 2)
  computed_var <- d$params$scale^2 * pi^2 / 3
  expect_equal(computed_var, 4, tolerance = 1e-6)
})

test_that("logistic negative mean", {
  d <- make_dist("logistic", mean = -5, var = 10)
  expect_equal(d$params$location, -5, tolerance = 1e-6)
})

test_that("logistic from two quantiles q1+q3", {
  d <- make_dist("logistic", q1 = 2, q3 = 8)
  expect_equal(d$q(0.25), 2, tolerance = 1e-6)
  expect_equal(d$q(0.75), 8, tolerance = 1e-6)
})

test_that("logistic from arbitrary quantiles", {
  d <- make_dist("logistic", quantiles = list(c(0.1, -5), c(0.9, 15)))
  expect_equal(d$q(0.1), -5, tolerance = 1e-4)
  expect_equal(d$q(0.9), 15, tolerance = 1e-4)
})

test_that("logistic from mean+q3", {
  d <- make_dist("logistic", mean = 5, q3 = 8)
  expect_equal(d$params$location, 5, tolerance = 1e-6)
  expect_equal(d$q(0.75), 8, tolerance = 1e-4)
})

test_that("logistic mean+median inconsistent raises", {
  expect_error(make_dist("logistic", mean = 5, median = 6), "must equal median")
})

test_that("logistic mean+median underdetermined raises", {
  expect_error(make_dist("logistic", mean = 5, median = 5), "additional constraint")
})

test_that("logistic from mode+iqr", {
  d <- make_dist("logistic", mode = 5, iqr = 4)
  expect_equal(d$q(0.75) - d$q(0.25), 4, tolerance = 1e-6)
  expect_equal(d$params$location, 5, tolerance = 1e-6)
})

test_that("logistic exists feasible", {
  expect_true(dist_exists("logistic", mean = 0, var = 3))
})

test_that("logistic exists any mean works", {
  expect_true(dist_exists("logistic", mean = -100, var = 0.01))
})

test_that("logistic exists zero var", {
  expect_false(dist_exists("logistic", mean = 0, var = 0))
})
