test_that("triangular from mean+var+mode", {
  d <- make_dist("triangular", mean = 5, var = 2, mode = 4)
  a <- d$params$a; b <- d$params$b; c <- d$params$c
  expect_equal((a + b + c) / 3, 5, tolerance = 1e-9)
  expect_equal((a^2 + b^2 + c^2 - a * b - a * c - b * c) / 18, 2, tolerance = 1e-9)
  expect_equal(c, 4)
})

test_that("triangular symmetric special case", {
  d <- make_dist("triangular", mean = 0, var = 1, mode = 0)
  expect_equal(d$params$a, -d$params$b, tolerance = 1e-9)
  expect_equal(d$params$c, 0)
})

test_that("triangular methods", {
  d <- make_dist("triangular", mean = 5, var = 2, mode = 4)
  expect_equal(d$p(d$params$a), 0, tolerance = 1e-12)
  expect_equal(d$p(d$params$b), 1, tolerance = 1e-12)
  expect_length(d$r(7), 7)
})
