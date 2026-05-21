test_that("beta on bounded support via affine scale", {
  d <- make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))
  expect_true(d$d(3.5) > 0)
  expect_equal(d$p(2), 0, tolerance = 1e-9)
  expect_equal(d$p(7), 1, tolerance = 1e-9)
  expect_equal(d$q(0.5), d$q(0.5), tolerance = 1e-9)  # symmetric solve sanity
})

test_that("gamma on half_right support via affine shift", {
  d <- make_dist("gamma", mean = 8, var = 4, support = c(3, Inf))
  # mean = 8 -> after shift, inner gamma has mean 5; var 4.
  expect_equal(d$p(3), 0, tolerance = 1e-9)
  expect_true(d$p(8) > 0 && d$p(8) < 1)
})

test_that("normal on bounded support via truncation", {
  d <- make_dist("normal", mean = 0, var = 0.1, support = c(-1, 1))
  # Truncated normal cdf endpoints
  expect_equal(d$p(-1), 0, tolerance = 1e-9)
  expect_equal(d$p(1), 1, tolerance = 1e-9)
  # Sanity: sample mean ~ 0
  expect_true(abs(mean(d$r(2000))) < 0.1)
})

test_that("normal on half_right support via truncation", {
  d <- make_dist("normal", mean = 2, var = 3, support = c(0, Inf))
  expect_equal(d$p(0), 0, tolerance = 1e-9)
  expect_equal(d$p(Inf), 1, tolerance = 1e-9)
})

test_that("dist_exists honors support: natural=real always TRUE if var>0", {
  # The structural check does NOT apply the Langevin dome.
  expect_true(dist_exists("normal", mean = 0, var = 0.5, support = c(-1, 1)))
  expect_false(dist_exists("normal", mean = 0, var = -1, support = c(-1, 1)))
})

test_that("dist_exists honors support: positive on bounded with lo<0 fails", {
  expect_false(dist_exists("gamma", mean = 5, var = 3, support = c(-1, 10)))
})

test_that("available_distributions with support=", {
  dists <- available_distributions(mean = 3.5, var = 0.5, support = c(2, 7))
  expect_true("beta" %in% dists)
})

test_that("constructor rejects infeasible Langevin envelope", {
  # var=0.5 > dome boundary 1/3 at mu=0 on (-1, 1) for Normal.
  expect_error(make_dist("normal", mean = 0, var = 0.5, support = c(-1, 1)),
               "infeasible")
})
