test_that("gumbel from mean+var", {
  d <- make_dist("gumbel", mean = 5, var = 3)
  beta <- sqrt(6 * 3) / pi
  gam <- -digamma(1)
  expect_equal(d$params$scale, beta, tolerance = 1e-12)
  expect_equal(d$params$location, 5 - beta * gam, tolerance = 1e-12)
})

test_that("gumbel from mean+mode", {
  d <- make_dist("gumbel", mean = 5, mode = 3)
  expect_equal(d$params$location, 3)
  # mean - mode = beta * EulerGamma
  expect_equal(d$params$scale, (5 - 3) / (-digamma(1)), tolerance = 1e-12)
})

test_that("gumbel d/p/q/r methods", {
  d <- make_dist("gumbel", mean = 0, var = pi^2 / 6)
  expect_equal(d$params$scale, 1, tolerance = 1e-12)
  expect_length(d$r(7), 7)
})

test_that("gumbel feasibility", {
  expect_true(dist_exists("gumbel", mean = 0, var = 1))
  expect_false(dist_exists("gumbel", mean = 0, var = -1))
})
