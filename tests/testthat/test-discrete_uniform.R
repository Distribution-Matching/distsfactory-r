test_that("discrete_uniform from mean+var", {
  d <- make_dist("discrete_uniform", mean = 5, var = 10)
  # var = ((b-a+1)^2-1)/12 = 10 -> (b-a+1)^2 = 121 -> b-a = 10
  expect_equal(d$params$max - d$params$min, 10)
  expect_equal((d$params$min + d$params$max) / 2, 5)
})

test_that("discrete_uniform methods", {
  d <- make_dist("discrete_uniform", mean = 0, var = 2)
  expect_equal(d$p(d$params$max), 1, tolerance = 1e-12)
  expect_equal(d$p(d$params$min - 1), 0, tolerance = 1e-12)
})

test_that("discrete_uniform alias", {
  d <- make_dist("dunif", mean = 5, var = 10)
  expect_equal(d$name, "discrete_uniform")
})
