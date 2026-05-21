test_that("fdist from mean+var", {
  d <- make_dist("fdist", mean = 1.5, var = 10)
  # df2 = 2*1.5/0.5 = 6; check df1 such that var formula holds
  d2 <- d$params$df2
  d1 <- d$params$df1
  expect_equal(d2, 6, tolerance = 1e-9)
  var_expected <- 2 * d2^2 * (d1 + d2 - 2) / (d1 * (d2 - 2)^2 * (d2 - 4))
  expect_equal(var_expected, 10, tolerance = 1e-6)
})

test_that("fdist feasibility", {
  expect_true(dist_exists("fdist", mean = 1.5, var = 10))
  expect_false(dist_exists("fdist", mean = 0.9, var = 5))  # mean must exceed 1
  expect_false(dist_exists("fdist", mean = 2.5, var = 5))  # df2 = 10/3 < 4 -> var undefined
})

test_that("fdist aliases", {
  expect_equal(make_dist("f", mean = 1.5, var = 10)$name, "fdist")
  expect_equal(make_dist("fisher", mean = 1.5, var = 10)$name, "fdist")
})

test_that("fdist d/p/q/r methods", {
  d <- make_dist("fdist", mean = 1.5, var = 10)
  expect_true(d$d(1) > 0)
  expect_equal(d$p(0), 0)
  expect_length(d$r(5), 5)
})
