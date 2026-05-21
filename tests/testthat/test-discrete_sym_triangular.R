test_that("discrete_sym_triangular from mean+var", {
  d <- make_dist("discrete_sym_triangular", mean = 5, var = 4)
  # var = n*(n+2)/6 = 4 -> n*(n+2) = 24 -> n = 4
  expect_equal(d$params$mu, 5)
  expect_equal(d$params$n, 4)
})

test_that("dst pmf normalises", {
  d <- make_dist("discrete_sym_triangular", mean = 5, var = 4)
  ks <- (d$params$mu - d$params$n):(d$params$mu + d$params$n)
  expect_equal(sum(d$d(ks)), 1, tolerance = 1e-12)
})

test_that("dst cdf endpoints", {
  d <- make_dist("discrete_sym_triangular", mean = 5, var = 4)
  expect_equal(d$p(d$params$mu - d$params$n - 1), 0)
  expect_equal(d$p(d$params$mu + d$params$n), 1)
})
