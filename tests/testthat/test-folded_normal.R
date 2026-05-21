test_that("folded_normal from mean+var", {
  d <- make_dist("folded_normal", mean = 2.5, var = 1.2)
  mu <- d$params$location; sigma <- d$params$scale
  m <- sigma * sqrt(2 / pi) * exp(-mu^2 / (2 * sigma^2)) +
       mu * (1 - 2 * pnorm(-mu / sigma))
  v <- mu^2 + sigma^2 - m^2
  expect_equal(m, 2.5, tolerance = 1e-6)
  expect_equal(v, 1.2, tolerance = 1e-6)
})

test_that("folded_normal feasibility", {
  expect_true(dist_exists("folded_normal", mean = 2.5, var = 1.2))
  # var/mean^2 capped at pi/2 - 1 ~ 0.571 (half-normal limit)
  expect_false(dist_exists("folded_normal", mean = 1, var = 10))
})

test_that("folded_normal methods", {
  d <- make_dist("folded_normal", mean = 2, var = 1)
  expect_equal(d$p(0), 0, tolerance = 1e-12)
  expect_length(d$r(5), 5)
})
