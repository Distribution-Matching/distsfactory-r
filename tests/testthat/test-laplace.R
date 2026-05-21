test_that("laplace from mean+var", {
  d <- make_dist("laplace", mean = 0, var = 2)
  expect_equal(d$params$location, 0)
  expect_equal(d$params$scale, 1, tolerance = 1e-12)  # var = 2*b^2 -> b=1
})

test_that("laplace from two quantiles", {
  d <- make_dist("laplace", q1 = -1, q3 = 1)
  expect_equal(d$q(0.25), -1, tolerance = 1e-9)
  expect_equal(d$q(0.75),  1, tolerance = 1e-9)
})

test_that("laplace feasibility", {
  expect_true(dist_exists("laplace", mean = 0, var = 1))
  expect_false(dist_exists("laplace", mean = 0, var = 0))
})

test_that("laplace d/p/q methods", {
  d <- make_dist("laplace", mean = 0, var = 2)
  expect_equal(d$p(0), 0.5, tolerance = 1e-12)
  expect_equal(d$q(0.5), 0, tolerance = 1e-12)
  expect_true(d$d(0) > 0)
  expect_length(d$r(5), 5)
})
