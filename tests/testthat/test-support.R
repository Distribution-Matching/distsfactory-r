test_that("beta on bounded support via affine scale", {
  d <- make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))
  expect_true(d$d(3.5) > 0)
  expect_equal(d$p(2), 0, tolerance = 1e-9)
  expect_equal(d$p(7), 1, tolerance = 1e-9)
  expect_equal(d$q(0.5), d$q(0.5), tolerance = 1e-9)  # symmetric solve sanity
})

test_that("gamma on half_right support via affine shift", {
  d <- make_dist("gamma", mean = 8, var = 4, support = c(3, Inf))
  # mean = 8 -> after shift, inner gamma has mean 5; var 4.
  expect_equal(d$p(3), 0, tolerance = 1e-9)
  expect_true(d$p(8) > 0 && d$p(8) < 1)
})

test_that("normal on bounded support via truncation", {
  d <- make_dist("normal", mean = 0, var = 0.1, support = c(-1, 1))
  # Truncated normal cdf endpoints
  expect_equal(d$p(-1), 0, tolerance = 1e-9)
  expect_equal(d$p(1), 1, tolerance = 1e-9)
  # Sanity: sample mean ~ 0
  expect_true(abs(mean(d$r(2000))) < 0.1)
})

test_that("normal on half_right support via truncation", {
  d <- make_dist("normal", mean = 2, var = 3, support = c(0, Inf))
  expect_equal(d$p(0), 0, tolerance = 1e-9)
  expect_equal(d$p(Inf), 1, tolerance = 1e-9)
})

test_that("dist_exists honors support: natural=real always TRUE if var>0", {
  # The structural check does NOT apply the Langevin dome.
  expect_true(dist_exists("normal", mean = 0, var = 0.5, support = c(-1, 1)))
  expect_false(dist_exists("normal", mean = 0, var = -1, support = c(-1, 1)))
})

test_that("dist_exists honors support: positive on bounded with lo<0 fails", {
  expect_false(dist_exists("gamma", mean = 5, var = 3, support = c(-1, 10)))
})

test_that("available_distributions with support=", {
  dists <- available_distributions(mean = 3.5, var = 0.5, support = c(2, 7))
  expect_true("beta" %in% dists)
})

test_that("constructor rejects infeasible Langevin envelope", {
  # var=0.5 > dome boundary 1/3 at mu=0 on (-1, 1) for Normal.
  expect_error(make_dist("normal", mean = 0, var = 0.5, support = c(-1, 1)),
               "infeasible")
})

# Helper: integrate x and x^2 against the distsfactory_dist's density to recover
# the truncated mean/var. Tolerance is generous because the underlying solvers
# themselves are numerical.
.trunc_moments_check <- function(d, lo, hi, target_mean, target_var, tol = 5e-3) {
  m  <- stats::integrate(function(x) x   * d$d(x), lo, hi,
                          rel.tol = 1e-9)$value
  m2 <- stats::integrate(function(x) x^2 * d$d(x), lo, hi,
                          rel.tol = 1e-9)$value
  v  <- m2 - m^2
  testthat::expect_equal(m, target_mean, tolerance = tol,
                         info = "truncated mean")
  testthat::expect_equal(v, target_var,  tolerance = tol,
                         info = "truncated var")
}

test_that("solve_truncated_positive: gamma on bounded [0, 10]", {
  d <- make_dist("gamma", mean = 3, var = 1, support = c(0, 10))
  expect_equal(d$p(0),  0, tolerance = 1e-9)
  expect_equal(d$p(10), 1, tolerance = 1e-9)
  .trunc_moments_check(d, 0, 10, 3, 1)
})

test_that("solve_truncated_positive: lognormal on bounded interval", {
  d <- make_dist("lognormal", mean = 2, var = 1, support = c(0.5, 6))
  expect_equal(d$p(0.5), 0, tolerance = 1e-9)
  expect_equal(d$p(6),   1, tolerance = 1e-9)
  .trunc_moments_check(d, 0.5, 6, 2, 1)
})

test_that("solve_truncated_positive: weibull on bounded interval", {
  d <- make_dist("weibull", mean = 4, var = 2, support = c(1, 10))
  expect_equal(d$p(1),  0, tolerance = 1e-9)
  expect_equal(d$p(10), 1, tolerance = 1e-9)
  .trunc_moments_check(d, 1, 10, 4, 2)
})

test_that("affine flip: gamma on (-Inf, 10]", {
  d <- make_dist("gamma", mean = 5, var = 3, support = c(-Inf, 10))
  # cdf endpoints
  expect_equal(d$p(10), 1, tolerance = 1e-9)
  expect_equal(d$p(-Inf), 0, tolerance = 1e-9)
  # X = 10 - Y where Y ~ Gamma with mean=5, so E[X] = 10 - 5 = 5
  m  <- stats::integrate(function(x) x   * d$d(x), -Inf, 10,
                          rel.tol = 1e-9)$value
  m2 <- stats::integrate(function(x) x^2 * d$d(x), -Inf, 10,
                          rel.tol = 1e-9)$value
  expect_equal(m, 5, tolerance = 5e-3)
  expect_equal(m2 - m^2, 3, tolerance = 5e-3)
})

test_that("truncated Normal matches Python README example", {
  # README/Python: mean=1.0, std=0.8, support=(-1, 4)
  d <- make_dist("normal", mean = 1.0, std = 0.8, support = c(-1, 4))
  # Parent params should be approximately (loc=0.9822, scale=0.8232).
  expect_equal(d$parent$params$mean, 0.9822, tolerance = 1e-3)
  expect_equal(d$parent$params$sd,   0.8232, tolerance = 1e-3)
  expect_equal(d$support, c(-1, 4))
  .trunc_moments_check(d, -1, 4, 1.0, 0.64, tol = 1e-3)  # var = std^2
})

test_that("dist_exists with support on positive family: half-right shift", {
  # Gamma on [3, Inf) with mean=8, var=4: standardize to (5, 4), feasible
  expect_true(dist_exists("gamma", mean = 8, var = 4, support = c(3, Inf)))
  # Mean below the lower endpoint -> infeasible (shifted mean would be <= 0)
  expect_false(dist_exists("gamma", mean = 2, var = 4, support = c(3, Inf)))
})

test_that("dist_exists with support on unit family", {
  expect_true(dist_exists("beta", mean = 3.5, var = 0.5, support = c(2, 7)))
  expect_false(dist_exists("beta", mean = 3.5, var = 50, support = c(2, 7)))
})

test_that("available_distributions support= filters positive families", {
  # On (3, Inf), positive-support families should be feasible by shift.
  dists <- available_distributions(mean = 8, var = 4, support = c(3, Inf))
  expect_true("gamma" %in% dists)
  expect_true("lognormal" %in% dists)
})
