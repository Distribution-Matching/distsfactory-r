test_that("erlang from mean+var (integer shape)", {
  d <- make_dist("erlang", mean = 4, var = 4)
  expect_equal(d$params$shape, 4)
  expect_equal(d$params$scale, 1)
})

test_that("erlang requires integer shape", {
  expect_error(make_dist("erlang", mean = 5, var = 3), "integer shape")
})

test_that("erlang feasibility", {
  expect_true(dist_exists("erlang", mean = 4, var = 4))
  expect_false(dist_exists("erlang", mean = 5, var = 3))
  expect_false(dist_exists("erlang", mean = -1, var = 1))
})
