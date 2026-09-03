# Coverage for marginal.plot(), which had none.
#
# Its output shape differs from partial.plot(): rather than a covariate grid,
# it returns one row per fitted subject in `$data`, plus a lowess of the same
# in `$smooth`, keyed by covariate and then by time point. It also reads the
# extrapolated `muhat` rather than the partial `mu`, so it needs its own
# exercise even though it shares a helper layer.

make.continuous <- function(n = 28) {
  d <- local({
    set.seed(3)
    simLong(
      n = n, n.time = 4, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  set.seed(7)
  list(
    fit = boostmtree(
      d$features, d$time, d$id, d$y,
      family = "continuous", M = 10, cv.flag = TRUE, verbose = FALSE
    ),
    d = d,
    n = n
  )
}

make.categorical <- function(family, y.fun) {
  d <- local({
    set.seed(8)
    simLong(
      n = 24, n.time = 3, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  # Build the response before seeding the fit. y.fun seeds itself, and R
  # forces arguments inside the call, so generating y in the argument list
  # would reset the RNG after set.seed() and leave the fit seed doing nothing.
  y <- y.fun(length(d$y))
  set.seed(7)
  boostmtree(
    d$features, d$time, d$id, y,
    family = family, M = 8, cv.flag = FALSE, verbose = FALSE
  )
}

test_that("marginal.plot returns per-subject data and a smooth, keyed by covariate then time", {
  o <- make.continuous()
  m <- marginal.plot(
    o$fit, x.var.names = c("x1", "x2"),
    time.points = c(0.5, 2), output = "data"
  )

  expect_s3_class(m, "marginal.plot.boostmtree")
  expect_named(m$data, c("x1", "x2"))
  expect_named(m$smooth, c("x1", "x2"))

  for (nm in c("x1", "x2")) {
    expect_length(m$data[[nm]], length(m$time.points))
    expect_length(m$smooth[[nm]], length(m$time.points))
    for (at.time in m$data[[nm]]) {
      expect_s3_class(at.time, "data.frame")
      expect_equal(ncol(at.time), 2L)
      expect_equal(nrow(at.time), o$n)
      expect_true(all(is.finite(as.matrix(at.time))))
    }
  }
})

test_that("the smooth tracks the data without reproducing it exactly", {
  o <- make.continuous()
  m <- marginal.plot(o$fit, x.var.names = "x1", time.points = 1, output = "data")
  raw <- m$data$x1[[1]]
  smoothed <- m$smooth$x1[[1]]

  expect_equal(nrow(smoothed), nrow(raw))
  expect_true(all(is.finite(as.matrix(smoothed))))
  expect_false(isTRUE(all.equal(raw[[2]], smoothed[[2]])))
  # a lowess is monotone in x, so its x column is the sorted covariate
  expect_false(is.unsorted(smoothed[[1]]))
})

test_that("subset restricts marginal.plot to the chosen subjects", {
  o <- make.continuous()
  keep <- o$fit$x$x1 > stats::median(o$fit$x$x1)
  m <- marginal.plot(
    o$fit, x.var.names = "x1", time.points = 1,
    subset = keep, output = "data"
  )
  expect_equal(nrow(m$data$x1[[1]]), sum(keep))
})

test_that("marginal.plot handles a binary response", {
  fit <- make.categorical("binary", function(n) {
    set.seed(9); rbinom(n, 1, 0.4)
  })
  m <- marginal.plot(fit, x.var.names = "x1", time.points = 1, output = "data")
  expect_true(all(is.finite(as.matrix(m$data$x1[[1]]))))
})

test_that("marginal.plot exposes one set per class for a nominal response", {
  fit <- make.categorical("nominal", function(n) {
    set.seed(9); factor(sample(c("a", "b", "c"), n, replace = TRUE))
  })
  m <- marginal.plot(fit, x.var.names = "x1", time.points = 1, output = "data")

  expect_gt(length(m$response.labels), 1L)
  expect_length(m$data, length(m$response.labels))
  for (per.response in m$data) {
    expect_true(all(is.finite(as.matrix(per.response$x1[[1]]))))
  }
})

test_that("an unmatched covariate name is rejected", {
  o <- make.continuous()
  expect_error(
    marginal.plot(o$fit, x.var.names = "not_a_variable", output = "data"),
    "does not match any fitted covariate"
  )
})

test_that("marginal plotting emits no warnings", {
  o <- make.continuous()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(marginal.plot(o$fit, x.var.names = "x1", time.points = 1))
})
