test_that("lognormal from mean+var basic", {
  d <- make_dist("lognormal", mean = 5, var = 3)
  s2 <- d$params$sdlog^2
  m_hat <- exp(d$params$meanlog + s2 / 2)
  v_hat <- (exp(s2) - 1) * exp(2 * d$params$meanlog + s2)
  expect_equal(m_hat, 5, tolerance = 1e-9)
  expect_equal(v_hat, 3, tolerance = 1e-9)
})

test_that("lognormal from mean+std", {
  d <- make_dist("lognormal", mean = 3, std = sqrt(2))
  s2 <- d$params$sdlog^2
  expect_equal((exp(s2) - 1) * exp(2 * d$params$meanlog + s2), 2, tolerance = 1e-9)
})

test_that("lognormal from two quantiles", {
  d <- make_dist("lognormal", quantiles = list(c(0.1, 1.0), c(0.9, 10.0)))
  expect_equal(d$q(0.1), 1.0, tolerance = 1e-6)
  expect_equal(d$q(0.9), 10.0, tolerance = 1e-6)
})

test_that("lognormal from mean+median", {
  d <- make_dist("lognormal", mean = 3, median = 2.5)
  s2 <- d$params$sdlog^2
  expect_equal(exp(d$params$meanlog + s2 / 2), 3, tolerance = 1e-6)
  expect_equal(d$q(0.5), 2.5, tolerance = 1e-6)
})

test_that("lognormal from mean+mode", {
  d <- make_dist("lognormal", mean = 5, mode = 2)
  s2 <- d$params$sdlog^2
  expect_equal(exp(d$params$meanlog + s2 / 2), 5, tolerance = 1e-9)
  expect_equal(exp(d$params$meanlog - s2), 2, tolerance = 1e-9)
})

test_that("lognormal exists predicate", {
  expect_true(dist_exists("lognormal", mean = 5, var = 3))
  expect_false(dist_exists("lognormal", mean = -1, var = 3))
  expect_false(dist_exists("lognormal", mean = 5, var = 0))
})

test_that("lognormal d/p/q/r methods", {
  d <- make_dist("lognormal", mean = 1, var = 0.25)
  expect_true(d$d(1) > 0)
  expect_equal(d$p(0), 0)
  expect_equal(d$p(Inf), 1)
  expect_length(d$r(20), 20)
})

test_that("lognormal alias lnorm", {
  d <- make_dist("lnorm", mean = 5, var = 3)
  expect_equal(d$name, "lognormal")
})
