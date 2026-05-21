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

# Roundtrip: build, recompute mean and var from the pmf directly, and verify
# the construction lands "reasonably close" to the requested (mean, var). The
# search is over integer triples in a +/-1 neighborhood of the continuous
# triangular solution (mirrors Python's implementation), so a perfect match is
# not expected — we just regression-protect the achieved accuracy.
.dtri_achieved <- function(params) {
  a <- params$a; b <- params$b; c <- params$c
  Z <- (b - a + 2) / 2
  ks <- a:b
  pmf <- ifelse(ks <= c,
                (ks - a + 1) / (c - a + 1) / Z,
                (b - ks + 1) / (b - c + 1) / Z)
  m <- sum(ks * pmf)
  list(mean = m, var = sum((ks - m)^2 * pmf), pmf_sum = sum(pmf))
}

test_that("dtri roundtrip: mode at peak with smallish var", {
  d <- make_dist("discrete_triangular", mean = 5, var = 2, mode = 5)
  ach <- .dtri_achieved(d$params)
  expect_equal(ach$pmf_sum, 1, tolerance = 1e-12)
  # Achieved (mean, var) within 1 absolute unit of target (loose bound — the
  # integer-triple search is necessarily imprecise; this only guards against
  # the search regressing wildly).
  expect_lt(abs(ach$mean - 5), 1.0)
  expect_lt(abs(ach$var - 2),  1.0)
})

test_that("dtri roundtrip: asymmetric (mode < mean)", {
  d <- make_dist("discrete_triangular", mean = 6, var = 3, mode = 4)
  expect_equal(d$params$c, 4)
  ach <- .dtri_achieved(d$params)
  expect_equal(ach$pmf_sum, 1, tolerance = 1e-12)
  expect_lt(abs(ach$mean - 6), 1.0)
  expect_lt(abs(ach$var - 3),  1.5)
})

test_that("dtri roundtrip: wide support", {
  d <- make_dist("discrete_triangular", mean = 10, var = 8, mode = 10)
  ach <- .dtri_achieved(d$params)
  expect_lt(abs(ach$mean - 10), 1.0)
  # Larger variance allows looser tolerance.
  expect_lt(abs(ach$var - 8),  2.0)
})

test_that("dtri q-method returns valid integer in support", {
  d <- make_dist("discrete_triangular", mean = 5, var = 2, mode = 5)
  qs <- d$q(c(0.1, 0.5, 0.9))
  expect_true(all(qs == round(qs)))
  expect_true(all(qs >= d$params$a & qs <= d$params$b))
})

test_that("dtri rejects non-positive variance", {
  expect_error(make_dist("discrete_triangular", mean = 5, var = 0, mode = 5),
               "variance must be positive")
  expect_error(make_dist("discrete_triangular", mean = 5, var = -1, mode = 5),
               "variance must be positive")
})
