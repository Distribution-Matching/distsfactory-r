test_that("binomial from mean+var", {
  d <- make_dist("binomial", mean = 5, var = 2.5)
  expect_equal(d$params$size, 10)
  expect_equal(d$params$prob, 0.5)
})

test_that("binomial requires integer n", {
  expect_error(make_dist("binomial", mean = 3.3, var = 2), "integer n")
})

test_that("binomial requires var < mean", {
  expect_error(make_dist("binomial", mean = 5, var = 5), "var < mean")
  expect_error(make_dist("binomial", mean = 5, var = 6), "var < mean")
})

test_that("binomial feasibility", {
  expect_true(dist_exists("binomial", mean = 5, var = 2.5))
  expect_false(dist_exists("binomial", mean = 5, var = 5))
})

test_that("binomial methods", {
  d <- make_dist("binomial", mean = 5, var = 2.5)
  expect_equal(d$d(5), dbinom(5, 10, 0.5))
  expect_equal(d$p(5), pbinom(5, 10, 0.5))
})
