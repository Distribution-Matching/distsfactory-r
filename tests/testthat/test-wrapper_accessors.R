# 0.2.0 wrapper-accessor surface. Every transform exposes $parent and
# $support uniformly; transform-specific metadata lives on a sibling slot:
#   - truncation: $support     = c(lo, hi)
#   - shift:      $shift       = a;  $support = c(a, Inf)
#   - flip:       $flip_point  = b;  $support = c(-Inf, b)
#   - scale:      $scale_loc   = a;  $scale_width = w;  $support = c(a, a+w)

# --- Truncated -----------------------------------------------------------

test_that("truncated normal: $parent is the unconstrained Normal", {
  d <- make_dist("normal", mean = 1, std = 0.8, support = c(-1, 4))
  expect_s3_class(d$parent, "distsfactory_dist")
  expect_equal(d$parent$name, "normal")
  expect_equal(d$parent$params$mean, 0.9822, tolerance = 1e-3)
  expect_equal(d$parent$params$sd,   0.8232, tolerance = 1e-3)
  expect_equal(d$support, c(-1, 4))
  expect_null(d$params)
})

test_that("truncated parent has its own moments (different from the truncated dist)", {
  d <- make_dist("normal", mean = 1, std = 0.8, support = c(-1, 4))
  expect_equal(d$mean(), 1.0, tolerance = 1e-6)
  expect_equal(d$parent$mean(), 0.9822, tolerance = 1e-3)
  expect_true(abs(d$mean() - d$parent$mean()) > 1e-3)  # genuinely different
})

# --- Shifted -------------------------------------------------------------

test_that("shifted gamma exposes $parent, $shift, $support", {
  d <- make_dist("gamma", mean = 8, var = 4, support = c(3, Inf))
  expect_s3_class(d$parent, "distsfactory_dist")
  expect_equal(d$parent$name, "gamma")
  expect_equal(d$parent$mean(), 5, tolerance = 1e-9)
  expect_equal(d$parent$var(), 4,  tolerance = 1e-9)
  expect_equal(d$shift, 3)
  expect_equal(d$support, c(3, Inf))
  expect_null(d$params)
})

test_that("shifted dist mean = parent mean + shift", {
  d <- make_dist("gamma", mean = 8, var = 4, support = c(3, Inf))
  expect_equal(d$mean(), d$parent$mean() + d$shift, tolerance = 1e-9)
  expect_equal(d$var(),  d$parent$var(),            tolerance = 1e-9)
})

# --- Flipped -------------------------------------------------------------

test_that("flipped gamma exposes $parent, $flip_point, $support", {
  d <- make_dist("gamma", mean = 5, var = 3, support = c(-Inf, 10))
  expect_s3_class(d$parent, "distsfactory_dist")
  expect_equal(d$parent$name, "gamma")
  expect_equal(d$flip_point, 10)
  expect_equal(d$support, c(-Inf, 10))
  expect_null(d$params)
})

test_that("flipped dist mean = flip_point - parent mean", {
  d <- make_dist("gamma", mean = 5, var = 3, support = c(-Inf, 10))
  expect_equal(d$mean(), d$flip_point - d$parent$mean(), tolerance = 1e-9)
  expect_equal(d$var(),  d$parent$var(),                  tolerance = 1e-9)
})

# --- Scaled --------------------------------------------------------------

test_that("scaled beta exposes $parent, $scale_loc, $scale_width, $support", {
  d <- make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))
  expect_s3_class(d$parent, "distsfactory_dist")
  expect_equal(d$parent$name, "beta")
  expect_equal(d$scale_loc, 2)
  expect_equal(d$scale_width, 5)
  expect_equal(d$support, c(2, 7))
  expect_null(d$params)
})

test_that("scaled dist mean/var transform correctly", {
  d <- make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))
  a <- d$scale_loc; w <- d$scale_width
  expect_equal(d$mean(), a + w * d$parent$mean(), tolerance = 1e-9)
  expect_equal(d$var(),  w^2 * d$parent$var(),    tolerance = 1e-9)
})

# --- print on wrappers ---------------------------------------------------

test_that("print on truncated dist mentions support and parent", {
  d <- make_dist("normal", mean = 1, std = 0.8, support = c(-1, 4))
  out <- capture.output(print(d))
  expect_match(paste(out, collapse = "\n"), "truncated_normal")
  expect_match(paste(out, collapse = "\n"), "truncated to \\[-1, 4\\]")
  expect_match(paste(out, collapse = "\n"), "parent: normal")
})

test_that("print on shifted dist mentions the shift and parent", {
  d <- make_dist("gamma", mean = 8, var = 4, support = c(3, Inf))
  out <- capture.output(print(d))
  expect_match(paste(out, collapse = "\n"), "shifted by 3")
  expect_match(paste(out, collapse = "\n"), "parent: gamma")
})

test_that("print on flipped dist mentions the flip point and parent", {
  d <- make_dist("gamma", mean = 5, var = 3, support = c(-Inf, 10))
  out <- capture.output(print(d))
  expect_match(paste(out, collapse = "\n"), "flipped about 10")
  expect_match(paste(out, collapse = "\n"), "parent: gamma")
})

test_that("print on scaled dist mentions the scaled support and parent", {
  d <- make_dist("beta", mean = 3.5, var = 0.5, support = c(2, 7))
  out <- capture.output(print(d))
  expect_match(paste(out, collapse = "\n"), "scaled to \\[2, 7\\]")
  expect_match(paste(out, collapse = "\n"), "parent: beta")
})

# --- $support always on wrappers, never on plain family dists -----------

test_that("plain family dist has no $support slot", {
  d <- make_dist("gamma", mean = 5, var = 3)
  expect_null(d$support)
  expect_null(d$parent)
})

test_that("every wrapped dist has a $support and a $parent", {
  cases <- list(
    make_dist("normal", mean = 1, std = 0.8, support = c(-1, 4)),
    make_dist("gamma",  mean = 8, var = 4,   support = c(3, Inf)),
    make_dist("gamma",  mean = 5, var = 3,   support = c(-Inf, 10)),
    make_dist("beta",   mean = 3.5, var = 0.5, support = c(2, 7))
  )
  for (d in cases) {
    expect_false(is.null(d$support), info = d$name)
    expect_false(is.null(d$parent),  info = d$name)
    expect_s3_class(d$parent, "distsfactory_dist")
  }
})
