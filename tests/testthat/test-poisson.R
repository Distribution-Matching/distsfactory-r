test_that("poisson from mean", {
  d <- make_dist("poisson", mean = 5)
  expect_equal(d$params$lambda, 5)
})

test_that("poisson rejects mean != var", {
  expect_error(make_dist("poisson", mean = 5, var = 3), "var == mean")
})

test_that("poisson feasibility", {
  expect_true(dist_exists("poisson", mean = 5, var = 5))
  expect_false(dist_exists("poisson", mean = 5, var = 3))
})

test_that("poisson methods", {
  d <- make_dist("poisson", mean = 3)
  expect_equal(d$d(2), dpois(2, 3))
  expect_length(d$r(7), 7)
})
