test_that("uniform from mean+var", {
  d <- make_dist("uniform", mean = 5, var = 3)
  half <- sqrt(9)
  expect_equal(d$params$min, 5 - half)
  expect_equal(d$params$max, 5 + half)
})

test_that("uniform from two quantiles", {
  d <- make_dist("uniform", q1 = 0, q3 = 4)
  # width = (4-0)/(0.75-0.25) = 8; a = 0 - 0.25*8 = -2; b = 6
  expect_equal(d$params$min, -2)
  expect_equal(d$params$max, 6)
})

test_that("uniform from mean+q1", {
  d <- make_dist("uniform", mean = 0, q1 = -2)
  # w = (-2 - 0)/(0.25 - 0.5) = 8; min = -4, max = 4
  expect_equal(d$params$min, -4)
  expect_equal(d$params$max,  4)
})

test_that("uniform feasibility", {
  expect_true(dist_exists("uniform", mean = 0, var = 1))
  expect_false(dist_exists("uniform", mean = 0, var = -1))
})

test_that("uniform d/p/q/r methods", {
  d <- make_dist("uniform", mean = 0.5, var = 1/12)  # standard uniform
  expect_equal(d$d(0.5), 1)
  expect_equal(d$p(0.5), 0.5)
  expect_equal(d$q(0.25), 0.25)
  expect_length(d$r(10), 10)
})

test_that("uniform alias unif", {
  d <- make_dist("unif", mean = 0, var = 1)
  expect_equal(d$name, "uniform")
})
