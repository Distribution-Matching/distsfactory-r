test_that("rayleigh from mean", {
  d <- make_dist("rayleigh", mean = 5)
  expect_equal(d$params$scale, 5 / sqrt(pi / 2), tolerance = 1e-12)
})

test_that("rayleigh from mean+var consistent", {
  d <- make_dist("rayleigh", mean = 5, var = 5^2 * (4 - pi) / pi)
  expect_equal(d$params$scale, 5 / sqrt(pi / 2), tolerance = 1e-9)
})

test_that("rayleigh from mode", {
  d <- make_dist("rayleigh", mode = 2)
  expect_equal(d$params$scale, 2)
})

test_that("rayleigh from median", {
  d <- make_dist("rayleigh", median = 2)
  expect_equal(d$params$scale, 2 / sqrt(2 * log(2)), tolerance = 1e-12)
})

test_that("rayleigh feasibility", {
  expect_true(dist_exists("rayleigh", mean = 5, var = 5^2 * (4 - pi) / pi))
  expect_false(dist_exists("rayleigh", mean = 5, var = 1))
})
