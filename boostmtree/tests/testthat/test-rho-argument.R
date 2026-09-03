# Characterization test for the `rho` argument.
#
# In the pre-2.0.0 layout the documented `rho` argument was silently ignored: a
# hidden `is.hidden.rho(user.option)` mechanism read `list(...)$rho`, but a named
# `rho =` call always binds to the formal, so the hidden lookup was always NULL
# and the formal was then overwritten with 0. The 2.0.0 refactor replaced that
# with `boostmtree.initialize.rho()`, which honours the formal.
#
# These tests pin that behaviour. They assert the three things that together
# mean "honoured" rather than merely "recorded": the value is held, it is
# validated, and it actually changes the fit.

make.data <- function() {
  set.seed(1)
  simLong(
    n = 60, n.time = 5, rho = 0.8, model = 2,
    family = "continuous", q = 0
  )$data.list
}

fit.with.rho <- function(d, rho) {
  set.seed(9)
  boostmtree(
    d$features, d$time, d$id, d$y,
    family = "continuous", M = 15, cv.flag = FALSE,
    verbose = FALSE, rho = rho
  )
}

test_that("a supplied rho is held fixed rather than re-estimated", {
  d <- make.data()
  for (value in c(0, 0.5, 0.9)) {
    fit <- fit.with.rho(d, value)
    expect_equal(unique(fit$rho), value)
  }
})

test_that("rho = NULL estimates rho rather than holding it fixed", {
  # Asserting "more than one distinct value across iterations" would be
  # over-specified: estimation can legitimately converge to a constant on some
  # data or platform. The property that actually matters is that NULL does not
  # behave like a fixed-rho run.
  d <- make.data()
  estimated <- fit.with.rho(d, NULL)$rho
  fixed <- fit.with.rho(d, 0)$rho

  expect_equal(length(estimated), length(fixed))
  expect_false(isTRUE(all.equal(estimated, fixed)))
})

test_that("rho outside (-1, 1) is rejected rather than silently coerced", {
  d <- make.data()
  for (bad in c(1, -1, 1.5)) {
    expect_error(fit.with.rho(d, bad), "strictly inside")
  }
})

test_that("rho materially changes the fit, not just the recorded value", {
  d <- make.data()
  independent <- unlist(fit.with.rho(d, 0)$mu)
  correlated <- unlist(fit.with.rho(d, 0.9)$mu)

  expect_equal(length(independent), length(correlated))
  expect_false(isTRUE(all.equal(independent, correlated)))
  expect_gt(max(abs(independent - correlated)), 1e-6)
})
