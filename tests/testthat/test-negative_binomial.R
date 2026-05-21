test_that("negative_binomial from mean+var", {
  d <- make_dist("negative_binomial", mean = 5, var = 8)
  expect_equal(d$params$prob, 5 / 8, tolerance = 1e-12)
  # r = mean * p / (1 - p) = 5 * 5/8 / (3/8) = 25/3
  expect_equal(d$params$size, 25 / 3, tolerance = 1e-12)
})

test_that("negative_binomial requires var > mean", {
  expect_error(make_dist("negative_binomial", mean = 5, var = 3), "var > mean")
})

test_that("negative_binomial feasibility", {
  expect_true(dist_exists("negative_binomial", mean = 5, var = 8))
  expect_false(dist_exists("negative_binomial", mean = 5, var = 3))
})
