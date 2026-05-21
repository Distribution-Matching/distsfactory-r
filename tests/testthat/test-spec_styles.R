# Coverage for spec-style combinations exercised by the README/Python sibling
# but not directly covered by the per-family tests. Each test verifies the
# documented input shape resolves to a distribution whose moments / quantiles
# match the targets.

# --- Moment alternatives ------------------------------------------------

test_that("second_moment is converted to var (gamma)", {
  d <- make_dist("gamma", mean = 5, second_moment = 28)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-9)
  expect_equal(d$params$shape / d$params$rate^2, 28 - 25, tolerance = 1e-9)
})

test_that("second_moment is converted to var (normal)", {
  d <- make_dist("normal", mean = 0, second_moment = 4)
  expect_equal(d$params$sd, 2)
})

test_that("std is converted to var (lognormal)", {
  d <- make_dist("lognormal", mean = 5, std = sqrt(3))
  m <- exp(d$params$meanlog + d$params$sdlog^2 / 2)
  s2 <- d$params$sdlog^2
  v <- (exp(s2) - 1) * exp(2 * d$params$meanlog + s2)
  expect_equal(m, 5, tolerance = 1e-9)
  expect_equal(v, 3, tolerance = 1e-9)
})

test_that("cv is converted to var (gamma)", {
  d <- make_dist("gamma", mean = 5, cv = 0.5)
  v_target <- (0.5 * 5)^2  # var = (cv*mean)^2 = 6.25
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-9)
  expect_equal(d$params$shape / d$params$rate^2, v_target, tolerance = 1e-9)
})

test_that("scv is converted to var (gamma)", {
  d <- make_dist("gamma", mean = 4, scv = 0.5)
  v_target <- 0.5 * 16
  expect_equal(d$params$shape / d$params$rate^2, v_target, tolerance = 1e-9)
})

# --- Quantile alternatives ----------------------------------------------

test_that("median + iqr resolves to two-quantile spec (normal)", {
  d <- make_dist("normal", median = 5, iqr = 2)
  # median = 5, iqr = 2 -> q1 = 4, q3 = 6 -> sigma = 1 / qnorm(0.75)
  expect_equal(d$params$mean, 5, tolerance = 1e-9)
  expect_equal(d$params$sd, 1 / qnorm(0.75), tolerance = 1e-9)
  expect_equal(d$q(0.25), 4, tolerance = 1e-9)
  expect_equal(d$q(0.75), 6, tolerance = 1e-9)
})

test_that("median + iqr (logistic)", {
  d <- make_dist("logistic", median = 0, iqr = 4)
  # Logistic: location=median, scale = iqr / (2*log(3))
  expect_equal(d$params$location, 0, tolerance = 1e-9)
  expect_equal(d$q(0.75) - d$q(0.25), 4, tolerance = 1e-9)
})

test_that("q1 alone with mean (mean+quantile spec)", {
  d <- make_dist("normal", mean = 0, q1 = -1)
  # sigma = (-1 - 0) / qnorm(0.25) = -1 / -0.6745 = 1.4826
  expect_equal(d$params$sd, -1 / qnorm(0.25), tolerance = 1e-9)
})

test_that("q3 alone with mean (mean+quantile spec)", {
  d <- make_dist("normal", mean = 0, q3 = 1)
  expect_equal(d$params$sd, 1 / qnorm(0.75), tolerance = 1e-9)
})

test_that("median alone (single-quantile, exponential)", {
  d <- make_dist("exponential", median = 2)
  # qexp(0.5, rate) = log(2)/rate = 2 -> rate = log(2)/2
  expect_equal(d$params$rate, log(2) / 2, tolerance = 1e-9)
})

test_that("single quantile (rayleigh median)", {
  d <- make_dist("rayleigh", median = 2)
  expect_equal(d$params$scale, 2 / sqrt(2 * log(2)), tolerance = 1e-9)
})

# --- Mode-augmented specs -----------------------------------------------

test_that("mode+var on Normal (symmetric: mode == mean)", {
  d <- make_dist("normal", mode = 3, var = 4)
  expect_equal(d$params$mean, 3)
  expect_equal(d$params$sd, 2)
})

test_that("mode+iqr on logistic", {
  d <- make_dist("logistic", mode = 5, iqr = 4)
  expect_equal(d$params$location, 5)  # logistic mode = location
  expect_equal(d$q(0.75) - d$q(0.25), 4, tolerance = 1e-9)
})

test_that("mode+quantile on gamma", {
  d <- make_dist("gamma", mode = 3, median = 4)
  expect_equal(d$q(0.5), 4, tolerance = 1e-4)
})

test_that("mode+iqr on gamma", {
  d <- make_dist("gamma", mode = 3, iqr = 4)
  expect_equal(d$q(0.75) - d$q(0.25), 4, tolerance = 1e-4)
})

test_that("mode alone (Rayleigh)", {
  d <- make_dist("rayleigh", mode = 2)
  expect_equal(d$params$scale, 2)
})

# --- 3-parameter triangular --------------------------------------------

test_that("triangular mean+var+mode roundtrip", {
  d <- make_dist("triangular", mean = 5, var = 2, mode = 4)
  expect_equal((d$params$a + d$params$b + d$params$c) / 3, 5, tolerance = 1e-9)
  expect_equal(d$params$c, 4)
})

# --- 1-parameter family from single moment -----------------------------

test_that("var alone (exponential)", {
  d <- make_dist("exponential", var = 9)
  expect_equal(d$params$rate, 1 / 3, tolerance = 1e-9)
})

test_that("var alone (geometric)", {
  d <- make_dist("geometric", var = 6)
  expect_equal(d$params$prob, 1 / 3, tolerance = 1e-9)
})

# --- Two-quantile on positive families ---------------------------------

test_that("two-quantile (lognormal)", {
  d <- make_dist("lognormal", quantiles = list(c(0.1, 1.0), c(0.9, 10.0)))
  expect_equal(d$q(0.1),  1.0, tolerance = 1e-6)
  expect_equal(d$q(0.9), 10.0, tolerance = 1e-6)
})

test_that("two-quantile (weibull)", {
  d <- make_dist("weibull", quantiles = list(c(0.1, 1.0), c(0.9, 5.0)))
  expect_equal(d$q(0.1), 1.0, tolerance = 1e-6)
  expect_equal(d$q(0.9), 5.0, tolerance = 1e-6)
})

test_that("two-quantile (pareto)", {
  d <- make_dist("pareto", quantiles = list(c(0.1, 1.0), c(0.9, 10.0)))
  expect_equal(d$q(0.1),  1.0, tolerance = 1e-6)
  expect_equal(d$q(0.9), 10.0, tolerance = 1e-6)
})

# --- Mean+quantile on positive families --------------------------------

test_that("mean+median on lognormal", {
  d <- make_dist("lognormal", mean = 3, median = 2.5)
  s2 <- d$params$sdlog^2
  expect_equal(exp(d$params$meanlog + s2 / 2), 3, tolerance = 1e-9)
  expect_equal(d$q(0.5), 2.5, tolerance = 1e-9)
})

test_that("mean+median on weibull", {
  d <- make_dist("weibull", mean = 3, median = 2.5)
  g1 <- gamma(1 + 1 / d$params$shape)
  expect_equal(d$params$scale * g1, 3, tolerance = 1e-6)
  expect_equal(d$q(0.5), 2.5, tolerance = 1e-6)
})

# --- Error paths --------------------------------------------------------

test_that("no spec raises", {
  expect_error(make_dist("gamma"), "at least one")
})

test_that("unsupported spec for family raises clearly", {
  # Gamma doesn't support a single-quantile spec.
  expect_error(make_dist("gamma", median = 4),
               "does not support specification type")
})

test_that("infeasible mean+var raises (exponential)", {
  expect_error(make_dist("exponential", mean = 2, var = 3),
               "must equal mean")
})
