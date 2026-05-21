test_that("partial_dist pinning gamma shape solves rate from mean", {
  spec <- partial_dist("gamma", shape = 3.0)
  d <- make_dist(spec, mean = 5.0)
  # mean = shape / rate -> rate = shape / mean = 3/5 = 0.6
  expect_equal(d$params$shape, 3.0)
  expect_equal(d$params$rate, 0.6, tolerance = 1e-9)
})

test_that("partial_dist pinning normal mean solves sd from var", {
  spec <- partial_dist("normal", mean = 1.0)
  d <- make_dist(spec, var = 4.0)
  expect_equal(d$params$mean, 1.0)
  expect_equal(d$params$sd, 2.0, tolerance = 1e-9)
})

test_that("partial_dist with no fixed params (2 free) solves mean+var", {
  spec <- partial_dist("gamma")  # nothing pinned
  d <- make_dist(spec, mean = 5, var = 3)
  expect_equal(d$params$shape / d$params$rate, 5, tolerance = 1e-6)
  expect_equal(d$params$shape / d$params$rate^2, 3, tolerance = 1e-6)
})

test_that("partial_dist normal alias / unknown param raises", {
  expect_error(partial_dist("gamma", nonsense = 1), "unknown parameter")
})

test_that("partial_dist with 1-param family + var-only solves", {
  spec <- partial_dist("exponential")  # rate is free, single param
  d <- make_dist(spec, mean = 5)
  expect_equal(d$params$rate, 0.2, tolerance = 1e-9)
})

test_that("partial_dist returns distsfactory_dist", {
  spec <- partial_dist("gamma", shape = 2)
  d <- make_dist(spec, mean = 4)
  expect_s3_class(d, "distsfactory_dist")
  expect_true(is.function(d$d))
  expect_true(is.function(d$p))
})

test_that("partial_dist print method", {
  spec <- partial_dist("gamma", shape = 2)
  expect_output(print(spec), "partial_dist\\(gamma")
})
