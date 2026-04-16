test_that("gamma from mean+var basic", {
  d <- make_dist("gamma", mean = 5, var = 3)
  expect_equal(d$params$shape * (1 / d$params$rate), 5, tolerance = 1e-6)
  # var = shape / rate^2
  expect_equal(d$params$shape / d$params$rate^2, 3, tolerance = 1e-6)
})

test_that("gamma from mean+var high variance", {
  d <- make_dist("gamma", mean = 10, var = 50)
  expect_equal(d$params$shape / d$params$rate, 10, tolerance = 1e-6)
  expect_equal(d$params$shape / d$params$rate^2, 50, tolerance = 1e-6)
})

test_that("gamma from mean+var low variance", {
  d <- make_dist("gamma", mean = 100, var = 1)
  expect_equal(d$params$shape / d$params$rate, 100, tolerance = 1e-6)
  expect_equal(d$params$shape / d$params$rate^2, 1, tolerance = 1e-6)
})

test_that("gamma from mean+std", {
  d <- make_dist("gamma", mean = 5, std = 2)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-6)
  expect_equal(d$params$shape / d$params$rate^2, 4, tolerance = 1e-6)
})

test_that("gamma from mean+cv", {
  d <- make_dist("gamma", mean = 5, cv = 0.5)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-6)
  expect_equal(d$params$shape / d$params$rate^2, 6.25, tolerance = 1e-6)
})

test_that("gamma from mean+scv", {
  d <- make_dist("gamma", mean = 4, scv = 0.5)
  expect_equal(d$params$shape / d$params$rate, 4, tolerance = 1e-6)
  expect_equal(d$params$shape / d$params$rate^2, 8, tolerance = 1e-6)
})

test_that("gamma from mean+mode basic", {
  d <- make_dist("gamma", mean = 5, mode = 3)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-6)
})

test_that("gamma from mean+mode close to zero", {
  d <- make_dist("gamma", mean = 5, mode = 0.5)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-6)
})

test_that("gamma mean+mode raises when mode > mean", {
  expect_error(make_dist("gamma", mean = 3, mode = 5), "must be greater than mode")
})

test_that("gamma from two quantiles q1+q3", {
  d <- make_dist("gamma", q1 = 2.0, q3 = 8.0)
  expect_equal(d$q(0.25), 2.0, tolerance = 1e-4)
  expect_equal(d$q(0.75), 8.0, tolerance = 1e-4)
})

test_that("gamma from arbitrary quantiles", {
  d <- make_dist("gamma", quantiles = list(c(0.1, 1.0), c(0.9, 10.0)))
  expect_equal(d$q(0.1), 1.0, tolerance = 1e-4)
  expect_equal(d$q(0.9), 10.0, tolerance = 1e-4)
})

test_that("gamma from mean+median", {
  d <- make_dist("gamma", mean = 5, median = 4.5)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-4)
  expect_equal(d$q(0.5), 4.5, tolerance = 1e-4)
})

test_that("gamma mode+median", {
  d <- make_dist("gamma", mode = 3, median = 4)
  expect_equal(d$q(0.5), 4.0, tolerance = 1e-4)
})

test_that("gamma mode+iqr", {
  d <- make_dist("gamma", mode = 3, iqr = 4)
  expect_equal(d$q(0.75) - d$q(0.25), 4.0, tolerance = 1e-4)
})

test_that("gamma exists feasible", {
  expect_true(dist_exists("gamma", mean = 5, var = 3))
})

test_that("gamma exists negative mean", {
  expect_false(dist_exists("gamma", mean = -1, var = 3))
})

test_that("gamma exists zero variance", {
  expect_false(dist_exists("gamma", mean = 5, var = 0))
})

test_that("gamma d/p/q/r methods work", {
  d <- make_dist("gamma", mean = 5, var = 3)
  expect_true(is.numeric(d$d(2)))
  expect_true(d$p(5) > 0 && d$p(5) < 1)
  expect_true(is.numeric(d$q(0.5)))
  expect_length(d$r(10), 10)
})
