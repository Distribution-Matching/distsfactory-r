test_that("frechet from mean+var", {
  d <- make_dist("frechet", mean = 5, var = 3)
  a <- d$params$shape; sg <- d$params$scale
  expect_equal(sg * gamma(1 - 1 / a), 5, tolerance = 1e-6)
  expect_equal(sg^2 * (gamma(1 - 2 / a) - gamma(1 - 1 / a)^2), 3, tolerance = 1e-6)
})

test_that("frechet feasibility", {
  expect_true(dist_exists("frechet", mean = 5, var = 3))
  expect_false(dist_exists("frechet", mean = -1, var = 3))
})

test_that("frechet methods", {
  d <- make_dist("frechet", mean = 5, var = 3)
  expect_true(d$d(1) >= 0)
  expect_length(d$r(5), 5)
})
