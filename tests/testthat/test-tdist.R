test_that("tdist from mean=0+var > 1", {
  d <- make_dist("tdist", mean = 0, var = 2)
  expect_equal(d$params$df, 4)  # var = 4/(4-2) = 2
})

test_that("tdist from var alone", {
  d <- make_dist("tdist", var = 5)
  expect_equal(d$params$df, 2.5)
})

test_that("tdist rejects nonzero mean", {
  expect_error(make_dist("tdist", mean = 1, var = 2), "mean 0")
})

test_that("tdist rejects var <= 1", {
  expect_error(make_dist("tdist", mean = 0, var = 1), "variance must exceed 1")
})

test_that("tdist feasibility", {
  expect_true(dist_exists("tdist", mean = 0, var = 2))
  expect_false(dist_exists("tdist", mean = 1, var = 2))
  expect_false(dist_exists("tdist", mean = 0, var = 0.5))
})

test_that("tdist alias 'student'/'t'", {
  d <- make_dist("student", var = 3)
  expect_equal(d$name, "tdist")
  d2 <- make_dist("t", var = 3)
  expect_equal(d2$name, "tdist")
})
