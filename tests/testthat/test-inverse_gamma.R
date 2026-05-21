test_that("inverse_gamma from mean+var", {
  d <- make_dist("inverse_gamma", mean = 5, var = 3)
  a <- d$params$shape; b <- d$params$scale
  expect_equal(b / (a - 1), 5, tolerance = 1e-9)
  expect_equal(b^2 / ((a - 1)^2 * (a - 2)), 3, tolerance = 1e-9)
})

test_that("inverse_gamma feasibility", {
  expect_true(dist_exists("inverse_gamma", mean = 5, var = 3))
  expect_false(dist_exists("inverse_gamma", mean = -1, var = 3))
})

test_that("inverse_gamma alias", {
  d <- make_dist("invgamma", mean = 5, var = 3)
  expect_equal(d$name, "inverse_gamma")
})

test_that("inverse_gamma methods", {
  d <- make_dist("inverse_gamma", mean = 5, var = 3)
  expect_true(d$d(1) >= 0)
  expect_length(d$r(5), 5)
})
