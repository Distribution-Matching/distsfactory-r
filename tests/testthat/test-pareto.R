test_that("pareto from mean+var", {
  d <- make_dist("pareto", mean = 5, var = 3)
  a <- d$params$shape; xm <- d$params$scale
  expect_equal(a * xm / (a - 1), 5, tolerance = 1e-9)
  expect_equal(xm^2 * a / ((a - 1)^2 * (a - 2)), 3, tolerance = 1e-9)
})

test_that("pareto from two quantiles", {
  d <- make_dist("pareto", quantiles = list(c(0.1, 1.0), c(0.9, 10.0)))
  expect_equal(d$q(0.1), 1.0, tolerance = 1e-6)
  expect_equal(d$q(0.9), 10.0, tolerance = 1e-6)
})

test_that("pareto feasibility", {
  expect_true(dist_exists("pareto", mean = 5, var = 3))
  expect_false(dist_exists("pareto", mean = -1, var = 3))
})

test_that("pareto methods", {
  d <- make_dist("pareto", mean = 5, var = 3)
  expect_true(d$d(d$params$scale + 1) > 0)
  expect_equal(d$p(d$params$scale), 0, tolerance = 1e-12)
  expect_length(d$r(5), 5)
})
