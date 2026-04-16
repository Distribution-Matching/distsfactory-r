test_that("beta from mean+var basic", {
  d <- make_dist("beta", mean = 0.4, var = 0.02)
  a <- d$params$shape1
  b <- d$params$shape2
  expect_equal(a / (a + b), 0.4, tolerance = 1e-6)
  computed_var <- a * b / ((a + b)^2 * (a + b + 1))
  expect_equal(computed_var, 0.02, tolerance = 1e-6)
})

test_that("beta from mean+var symmetric", {
  d <- make_dist("beta", mean = 0.5, var = 0.05)
  a <- d$params$shape1
  b <- d$params$shape2
  expect_equal(a / (a + b), 0.5, tolerance = 1e-6)
})

test_that("beta from mean+var skewed right", {
  d <- make_dist("beta", mean = 0.8, var = 0.01)
  a <- d$params$shape1
  b <- d$params$shape2
  expect_equal(a / (a + b), 0.8, tolerance = 1e-6)
})

test_that("beta from mean+std", {
  d <- make_dist("beta", mean = 0.5, std = 0.1)
  a <- d$params$shape1
  b <- d$params$shape2
  computed_var <- a * b / ((a + b)^2 * (a + b + 1))
  expect_equal(computed_var, 0.01, tolerance = 1e-6)
})

test_that("beta var too large raises", {
  expect_error(make_dist("beta", mean = 0.5, var = 0.3), "variance too large")
})

test_that("beta from mean+mode basic", {
  d <- make_dist("beta", mean = 0.4, mode = 0.35)
  a <- d$params$shape1
  b <- d$params$shape2
  expect_equal(a / (a + b), 0.4, tolerance = 1e-4)
})

test_that("beta mean+mode symmetric raises", {
  expect_error(make_dist("beta", mean = 0.5, mode = 0.5))
})

test_that("beta from two quantiles q1+q3", {
  d <- make_dist("beta", q1 = 0.2, q3 = 0.6)
  expect_equal(d$q(0.25), 0.2, tolerance = 1e-3)
  expect_equal(d$q(0.75), 0.6, tolerance = 1e-3)
})

test_that("beta from arbitrary quantiles", {
  d <- make_dist("beta", quantiles = list(c(0.1, 0.15), c(0.9, 0.85)))
  expect_equal(d$q(0.1), 0.15, tolerance = 1e-3)
  expect_equal(d$q(0.9), 0.85, tolerance = 1e-3)
})

test_that("beta from mean+median", {
  d <- make_dist("beta", mean = 0.4, median = 0.38)
  a <- d$params$shape1
  b <- d$params$shape2
  expect_equal(a / (a + b), 0.4, tolerance = 1e-3)
  expect_equal(d$q(0.5), 0.38, tolerance = 1e-3)
})

test_that("beta exists feasible", {
  expect_true(dist_exists("beta", mean = 0.5, var = 0.1))
})

test_that("beta exists var too large", {
  expect_false(dist_exists("beta", mean = 0.5, var = 0.3))
})

test_that("beta exists mean out of range", {
  expect_false(dist_exists("beta", mean = 1.5, var = 0.1))
  expect_false(dist_exists("beta", mean = -0.1, var = 0.1))
})

test_that("beta boundary var", {
  expect_false(dist_exists("beta", mean = 0.5, var = 0.25))
  expect_true(dist_exists("beta", mean = 0.5, var = 0.24))
})
