test_that("unknown distribution raises", {
  expect_error(make_dist("pareto", mean = 5), "Unknown distribution")
})

test_that("no spec raises", {
  expect_error(make_dist("gamma"), "at least one")
})

test_that("unsupported spec raises", {
  expect_error(make_dist("gamma", mean = 5), "does not support")
})

test_that("available_distributions mean+var", {
  dists <- available_distributions(mean = 5, var = 3)
  expect_true("gamma" %in% dists)
  expect_true("logistic" %in% dists)
})

test_that("available_distributions exponential only when var matches", {
  dists_match <- available_distributions(mean = 5, var = 25)
  expect_true("exponential" %in% dists_match)

  dists_no <- available_distributions(mean = 5, var = 3)
  expect_false("exponential" %in% dists_no)
})

test_that("available_distributions beta only in unit", {
  dists_unit <- available_distributions(mean = 0.5, var = 0.05)
  expect_true("beta" %in% dists_unit)

  dists_nope <- available_distributions(mean = 5, var = 3)
  expect_false("beta" %in% dists_nope)
})

test_that("dist_exists returns logical", {
  expect_true(dist_exists("gamma", mean = 5, var = 3))
  expect_type(dist_exists("gamma", mean = 5, var = 3), "logical")
})

test_that("dist_exists returns FALSE", {
  expect_false(dist_exists("exponential", mean = 5, var = 3))
})

test_that("case insensitive", {
  d <- make_dist("Gamma", mean = 5, var = 3)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-4)
})

test_that("mixed case", {
  d <- make_dist("BETA", mean = 0.5, var = 0.05)
  expect_equal(d$params$shape1 / (d$params$shape1 + d$params$shape2), 0.5, tolerance = 1e-4)
})

test_that("print method works", {
  d <- make_dist("gamma", mean = 5, var = 3)
  expect_output(print(d), "distsfactory: gamma")
})
