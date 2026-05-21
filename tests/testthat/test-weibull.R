test_that("weibull from mean+var", {
  d <- make_dist("weibull", mean = 5, var = 3)
  g1 <- gamma(1 + 1 / d$params$shape)
  g2 <- gamma(1 + 2 / d$params$shape)
  expect_equal(d$params$scale * g1, 5, tolerance = 1e-6)
  expect_equal(d$params$scale^2 * (g2 - g1^2), 3, tolerance = 1e-6)
})

test_that("weibull from two quantiles", {
  d <- make_dist("weibull", quantiles = list(c(0.1, 1.0), c(0.9, 5.0)))
  expect_equal(d$q(0.1), 1.0, tolerance = 1e-6)
  expect_equal(d$q(0.9), 5.0, tolerance = 1e-6)
})

test_that("weibull from mean+quantile", {
  d <- make_dist("weibull", mean = 3, median = 2.5)
  g1 <- gamma(1 + 1 / d$params$shape)
  expect_equal(d$params$scale * g1, 3, tolerance = 1e-6)
  expect_equal(d$q(0.5), 2.5, tolerance = 1e-6)
})

test_that("weibull feasibility", {
  expect_true(dist_exists("weibull", mean = 5, var = 3))
  expect_false(dist_exists("weibull", mean = -1, var = 3))
  expect_false(dist_exists("weibull", mean = 5, var = 0))
})

test_that("weibull d/p/q/r methods", {
  d <- make_dist("weibull", mean = 2, var = 1)
  expect_true(d$d(1) > 0)
  expect_equal(d$p(0), 0)
  expect_length(d$r(15), 15)
})
