test_that("chisq from mean only", {
  d <- make_dist("chisq", mean = 4)
  expect_equal(d$params$df, 4)
})

test_that("chisq from mean+var consistent", {
  d <- make_dist("chisq", mean = 4, var = 8)
  expect_equal(d$params$df, 4)
})

test_that("chisq rejects var != 2*mean", {
  expect_error(make_dist("chisq", mean = 4, var = 3), "var = 2.mean")
})

test_that("chisq feasibility", {
  expect_true(dist_exists("chisq", mean = 4, var = 8))
  expect_false(dist_exists("chisq", mean = 4, var = 9))
  expect_false(dist_exists("chisq", mean = 0, var = 0))
})

test_that("chisq d/p/q/r methods", {
  d <- make_dist("chisq", mean = 3)
  expect_true(d$d(1) > 0)
  expect_equal(d$p(0), 0)
  expect_length(d$r(5), 5)
})

test_that("chisq aliases", {
  expect_equal(make_dist("chi_sq", mean = 2)$name, "chisq")
  expect_equal(make_dist("chisquare", mean = 2)$name, "chisq")
})
