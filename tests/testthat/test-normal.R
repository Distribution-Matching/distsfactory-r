test_that("normal from mean+var basic", {
  d <- make_dist("normal", mean = 0, var = 1)
  expect_equal(d$params$mean, 0)
  expect_equal(d$params$sd, 1)
})

test_that("normal from mean+std", {
  d <- make_dist("normal", mean = 3, std = 2)
  expect_equal(d$params$mean, 3)
  expect_equal(d$params$sd, 2)
})

test_that("normal from mean+cv", {
  d <- make_dist("normal", mean = 4, cv = 0.5)
  expect_equal(d$params$sd, 2, tolerance = 1e-9)
})

test_that("normal from two quantiles", {
  d <- make_dist("normal", q1 = -1, q3 = 1)
  expect_equal(d$q(0.25), -1, tolerance = 1e-9)
  expect_equal(d$q(0.75),  1, tolerance = 1e-9)
})

test_that("normal mean+median (p=0.5) underdetermines sigma", {
  expect_error(make_dist("normal", mean = 5, median = 5), "quantile at p=0.5")
})

test_that("normal from mean+q1", {
  d <- make_dist("normal", mean = 0, q1 = -1)
  # sigma = (-1 - 0) / qnorm(0.25)
  sig <- -1 / qnorm(0.25)
  expect_equal(d$params$sd, sig, tolerance = 1e-9)
})

test_that("normal from mode+var (mode == mean)", {
  d <- make_dist("normal", mode = 3, var = 4)
  expect_equal(d$params$mean, 3)
  expect_equal(d$params$sd, 2)
})

test_that("normal d/p/q/r methods", {
  d <- make_dist("normal", mean = 0, var = 1)
  expect_equal(d$d(0), dnorm(0), tolerance = 1e-12)
  expect_equal(d$p(1), pnorm(1), tolerance = 1e-12)
  expect_equal(d$q(0.5), 0)
  expect_length(d$r(7), 7)
})

test_that("normal exists for any var > 0", {
  expect_true(dist_exists("normal", mean = 0, var = 1))
  expect_true(dist_exists("normal", mean = -100, var = 0.01))
  expect_false(dist_exists("normal", mean = 0, var = 0))
  expect_false(dist_exists("normal", mean = 0, var = -1))
})

test_that("normal alias 'norm' works", {
  d <- make_dist("norm", mean = 0, var = 1)
  expect_equal(d$params$sd, 1)
})
