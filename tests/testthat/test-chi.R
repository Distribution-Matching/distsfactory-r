test_that("chi from mean+var", {
  # df=4 has known mean = sqrt(2)*Gamma(2.5)/Gamma(2)
  m_target <- sqrt(2) * exp(lgamma(2.5) - lgamma(2))
  v_target <- 4 - m_target^2
  d <- make_dist("chi", mean = m_target, var = v_target)
  expect_equal(d$params$df, 4, tolerance = 1e-5)
})

test_that("chi from mean alone", {
  # E[Chi(3)] = sqrt(2)*Gamma(2)/Gamma(1.5) = sqrt(2) / (sqrt(pi)/2) = 2*sqrt(2/pi)
  target <- sqrt(2) * exp(lgamma(2) - lgamma(1.5))
  d <- make_dist("chi", mean = target)
  expect_equal(d$params$df, 3, tolerance = 1e-5)
})

test_that("chi feasibility", {
  m_target <- sqrt(2) * exp(lgamma(2.5) - lgamma(2))
  v_target <- 4 - m_target^2
  expect_true(dist_exists("chi", mean = m_target, var = v_target))
  expect_false(dist_exists("chi", mean = 1, var = 1))
})
