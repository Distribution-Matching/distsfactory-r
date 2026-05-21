test_that("sym_triangular from mean+var", {
  d <- make_dist("sym_triangular", mean = 2, var = 6)
  expect_equal(d$params$location, 2)
  expect_equal(d$params$scale, 6)  # a = sqrt(6*6) = 6
})

test_that("sym_triangular cdf at endpoints", {
  d <- make_dist("sym_triangular", mean = 0, var = 1)
  a <- d$params$scale
  expect_equal(d$p(-a), 0, tolerance = 1e-12)
  expect_equal(d$p(a), 1, tolerance = 1e-12)
  expect_equal(d$p(0), 0.5, tolerance = 1e-12)
})

test_that("sym_triangular quantile inverse", {
  d <- make_dist("sym_triangular", mean = 0, var = 1)
  expect_equal(d$p(d$q(0.3)), 0.3, tolerance = 1e-9)
})
