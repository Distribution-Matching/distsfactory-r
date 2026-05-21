test_that("discrete_triangular from mean+var+mode", {
  d <- make_dist("discrete_triangular", mean = 5, var = 2, mode = 5)
  expect_equal(d$params$c, 5)
  expect_true(d$params$a <= 5 && 5 <= d$params$b)
})

test_that("dtri pmf normalises", {
  d <- make_dist("discrete_triangular", mean = 5, var = 2, mode = 5)
  ks <- d$params$a:d$params$b
  expect_equal(sum(d$d(ks)), 1, tolerance = 1e-10)
})

test_that("dtri cdf endpoints", {
  d <- make_dist("discrete_triangular", mean = 5, var = 2, mode = 5)
  expect_equal(d$p(d$params$a - 1), 0)
  expect_equal(d$p(d$params$b), 1)
})
