test_that("cauchy from two quantiles q1 + q3", {
  d <- make_dist("cauchy", q1 = -1, q3 = 1)
  # Symmetric -> location = 0, scale = iqr/2 = 1
  expect_equal(d$params$location, 0, tolerance = 1e-9)
  expect_equal(d$params$scale, 1, tolerance = 1e-9)
  expect_equal(d$q(0.25), -1, tolerance = 1e-9)
  expect_equal(d$q(0.75),  1, tolerance = 1e-9)
})

test_that("cauchy from arbitrary two quantiles", {
  d <- make_dist("cauchy", quantiles = list(c(0.1, -3.0), c(0.9, 3.0)))
  expect_equal(d$q(0.1), -3.0, tolerance = 1e-9)
  expect_equal(d$q(0.9),  3.0, tolerance = 1e-9)
})

test_that("cauchy from median + iqr", {
  d <- make_dist("cauchy", median = 5, iqr = 4)
  expect_equal(d$params$location, 5, tolerance = 1e-9)
  expect_equal(d$params$scale, 2, tolerance = 1e-9)
})

test_that("cauchy rejects mean+var (no finite moments)", {
  expect_error(make_dist("cauchy", mean = 0, var = 1),
               "no finite mean or variance")
})

test_that("cauchy rejects mean-only", {
  expect_error(make_dist("cauchy", mean = 5),
               "no finite mean")
})

test_that("cauchy dist_exists is always FALSE for moment specs", {
  expect_false(dist_exists("cauchy", mean = 0, var = 1))
  expect_false(dist_exists("cauchy", mean = 5, var = 100))
})

test_that("cauchy d/p/q/r methods (matching base R)", {
  d <- make_dist("cauchy", q1 = -1, q3 = 1)
  expect_equal(d$d(0), dcauchy(0), tolerance = 1e-12)
  expect_equal(d$p(0), 0.5, tolerance = 1e-12)
  expect_equal(d$q(0.5), 0, tolerance = 1e-9)
  expect_length(d$r(7), 7)
})

test_that("cauchy d honors log; p honors lower.tail and log.p", {
  d <- make_dist("cauchy", median = 0, iqr = 2)
  expect_equal(d$d(0, log = TRUE), log(d$d(0)), tolerance = 1e-12)
  expect_equal(d$p(0, lower.tail = FALSE), 0.5, tolerance = 1e-12)
  expect_equal(d$p(0, log.p = TRUE), log(0.5), tolerance = 1e-12)
})

test_that("cauchy via partial_dist (pin location, solve from quantile)", {
  # partial_dist doesn't support quantile specs in 0.1.0, but pinning all
  # canonical params should round-trip.
  spec <- partial_dist("cauchy", location = 0, scale = 2)
  d <- make_dist(spec)
  expect_equal(d$params$location, 0)
  expect_equal(d$params$scale, 2)
})

test_that("cauchy alias 'lorentz'-style not exposed; only canonical works", {
  # Sanity: the canonical name resolves.
  d <- make_dist("Cauchy", q1 = -1, q3 = 1)   # case-insensitive
  expect_equal(d$name, "cauchy")
})
