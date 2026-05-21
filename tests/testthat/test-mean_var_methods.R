# `d$mean()`, `d$var()`, and `mean(d)` mirror the Julia (`mean(d)`, `var(d)`)
# and Python (`d.mean()`, `d.var()`) sibling-package conventions.

test_that("d$mean() and d$var() on a plain gamma match closed-form", {
  d <- make_dist("gamma", mean = 5, var = 3)
  expect_equal(d$mean(), 5, tolerance = 1e-9)
  expect_equal(d$var(), 3, tolerance = 1e-9)
})

test_that("mean(d) S3 method dispatches", {
  d <- make_dist("normal", mean = 2, var = 4)
  expect_equal(mean(d), 2)
})

test_that("d$mean() on a truncated normal returns the truncated mean", {
  # README example: mean=1, std=0.8 on [-1, 4]. Parent (mu, sigma) was solved
  # so the *truncated* mean is 1.0, NOT mu ≈ 0.982.
  d <- make_dist("normal", mean = 1.0, std = 0.8, support = c(-1, 4))
  expect_equal(d$mean(), 1.0, tolerance = 1e-6)
  expect_equal(d$var(), 0.64, tolerance = 1e-6)  # std^2
  # The parent Normal's parameters are reachable via $parent.
  expect_equal(d$parent$params$mean, 0.9822, tolerance = 1e-3)
  expect_equal(d$parent$params$sd,   0.8232, tolerance = 1e-3)
  expect_equal(d$support, c(-1, 4))
  # The wrapper itself has no $params (would be a category error).
  expect_null(d$params)
})

test_that("d$mean() on truncated half-right normal", {
  d <- make_dist("normal", mean = 2, var = 3, support = c(0, Inf))
  expect_equal(d$mean(), 2, tolerance = 1e-6)
  expect_equal(d$var(), 3, tolerance = 1e-6)
})

test_that("d$mean() on affine-shifted gamma", {
  # Gamma on [3, Inf) with target mean=8, var=4.
  d <- make_dist("gamma", mean = 8, var = 4, support = c(3, Inf))
  expect_equal(d$mean(), 8, tolerance = 1e-9)
  expect_equal(d$var(), 4, tolerance = 1e-9)
})

test_that("d$mean() on affine-flipped gamma", {
  d <- make_dist("gamma", mean = 5, var = 3, support = c(-Inf, 10))
  expect_equal(d$mean(), 5, tolerance = 1e-9)
  expect_equal(d$var(), 3, tolerance = 1e-9)
})

test_that("d$mean() on affine-scaled beta", {
  d <- make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))
  expect_equal(d$mean(), 3.5, tolerance = 1e-9)
  expect_equal(d$var(), 0.5, tolerance = 1e-9)
})

test_that("d$mean() returns NaN for Cauchy (no finite moments)", {
  d <- make_dist("cauchy", q1 = -1, q3 = 1)
  expect_true(is.nan(d$mean()))
  expect_true(is.nan(d$var()))
})

test_that("d$mean() and d$var() across a representative sample of families", {
  cases <- list(
    list(name = "lognormal",   args = list(mean = 5, var = 3)),
    list(name = "weibull",     args = list(mean = 5, var = 3)),
    list(name = "laplace",     args = list(mean = 0, var = 2)),
    list(name = "rayleigh",    args = list(mean = 5)),
    list(name = "binomial",    args = list(mean = 5, var = 2.5)),
    list(name = "poisson",     args = list(mean = 5)),
    list(name = "geometric",   args = list(mean = 2)),
    list(name = "triangular",  args = list(mean = 5, var = 2, mode = 4))
  )
  for (case in cases) {
    d <- do.call(make_dist, c(list(case$name), case$args))
    if (!is.null(case$args$mean)) {
      expect_equal(d$mean(), case$args$mean, tolerance = 1e-6,
                   info = sprintf("%s mean", case$name))
    }
    if (!is.null(case$args$var)) {
      expect_equal(d$var(), case$args$var, tolerance = 1e-6,
                   info = sprintf("%s var", case$name))
    }
  }
})
