test_that("geometric from mean", {
  d <- make_dist("geometric", mean = 2)
  expect_equal(d$params$prob, 1 / 3, tolerance = 1e-12)
})

test_that("geometric from mean+var consistent", {
  d <- make_dist("geometric", mean = 2, var = 2 * 3)
  expect_equal(d$params$prob, 1 / 3, tolerance = 1e-9)
})

test_that("geometric from var alone", {
  d <- make_dist("geometric", var = 6)
  expect_equal(d$params$prob, 1 / 3, tolerance = 1e-12)
})

test_that("geometric from quantile", {
  d <- make_dist("geometric", median = 1)
  # P(X <= 1) = 0.5 -> prob = 1 - 0.5^(1/2)
  expect_equal(d$params$prob, 1 - sqrt(0.5), tolerance = 1e-12)
})

test_that("geometric feasibility", {
  expect_true(dist_exists("geometric", mean = 2, var = 6))
  expect_false(dist_exists("geometric", mean = 2, var = 4))
})
