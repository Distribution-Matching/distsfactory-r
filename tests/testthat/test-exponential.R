test_that("exponential from mean basic", {
  d <- make_dist("exponential", mean = 3)
  expect_equal(1 / d$params$rate, 3, tolerance = 1e-6)
})

test_that("exponential from small mean", {
  d <- make_dist("exponential", mean = 0.1)
  expect_equal(1 / d$params$rate, 0.1, tolerance = 1e-6)
})

test_that("exponential from large mean", {
  d <- make_dist("exponential", mean = 1000)
  expect_equal(1 / d$params$rate, 1000, tolerance = 1e-6)
})

test_that("exponential from mean+var consistent", {
  d <- make_dist("exponential", mean = 2.5, var = 6.25)
  expect_equal(1 / d$params$rate, 2.5, tolerance = 1e-6)
})

test_that("exponential from mean+var inconsistent raises", {
  expect_error(make_dist("exponential", mean = 2.5, var = 1.5), "must equal mean")
})

test_that("exponential from var", {
  d <- make_dist("exponential", var = 9)
  expect_equal(1 / d$params$rate, 3, tolerance = 1e-6)
})

test_that("exponential from median", {
  d <- make_dist("exponential", median = 2.0)
  expect_equal(d$q(0.5), 2.0, tolerance = 1e-6)
})

test_that("exponential from q1", {
  d <- make_dist("exponential", q1 = 1.0)
  expect_equal(d$q(0.25), 1.0, tolerance = 1e-6)
})

test_that("exponential from q3", {
  d <- make_dist("exponential", q3 = 5.0)
  expect_equal(d$q(0.75), 5.0, tolerance = 1e-6)
})

test_that("exponential exists feasible", {
  expect_true(dist_exists("exponential", mean = 2.5, var = 6.25))
})

test_that("exponential exists infeasible var", {
  expect_false(dist_exists("exponential", mean = 2.5, var = 1.5))
})

test_that("exponential exists negative mean", {
  expect_false(dist_exists("exponential", mean = -1, var = 1))
})

test_that("exp alias works", {
  d <- make_dist("exp", mean = 3)
  expect_equal(1 / d$params$rate, 3, tolerance = 1e-6)
})

test_that("expon alias works", {
  d <- make_dist("expon", mean = 3)
  expect_equal(1 / d$params$rate, 3, tolerance = 1e-6)
})
