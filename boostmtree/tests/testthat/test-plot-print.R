# Coverage for the display layer: plot.boostmtree(), print.boostmtree(),
# boostmtree.news(), and the shared helpers in boostmtree_display.R. All were
# at 0%.
#
# Two v2.0.1 bugs (a truncated vimp label, and a scalar Mopt indexed as if it
# were a vector) both reached a release because nothing asserted the absence of
# warnings around a plotting call. These tests assert it.

fit.continuous <- function(cv.flag = TRUE, n = 26) {
  d <- local({
    set.seed(2)
    simLong(
      n = n, n.time = 4, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  set.seed(7)
  list(
    fit = boostmtree(
      d$features, d$time, d$id, d$y,
      family = "continuous", M = 10, cv.flag = cv.flag, verbose = FALSE
    ),
    d = d
  )
}

# Run an expression with a null graphics device, closing it even on failure.
# Deliberately not using withr::defer here: withr is not a declared dependency
# of this package, and reaching for it would make R CMD check flag an
# undeclared import.
on.null.device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

test_that("plot(output = 'data') returns finite per-iteration series", {
  o <- fit.continuous()
  p <- plot(o$fit, output = "data")

  expect_type(p, "list")
  expect_gt(length(p), 0L)
  expect_true(all(vapply(p, function(e) all(is.finite(unlist(e[sapply(e, is.numeric)]))), logical(1))))
})

test_that("use.rmse toggles the plotted quantity", {
  o <- fit.continuous()
  with.rmse <- plot(o$fit, output = "data", use.rmse = TRUE)
  without <- plot(o$fit, output = "data", use.rmse = FALSE)

  expect_equal(length(with.rmse), length(without))
  expect_false(isTRUE(all.equal(with.rmse, without)))
})

test_that("plotting a fit emits no warnings, with and without cross-validation", {
  on.null.device({
    expect_no_warning(plot(fit.continuous(cv.flag = TRUE)$fit))
    expect_no_warning(plot(fit.continuous(cv.flag = FALSE)$fit))
  })
})

test_that("plotting a predict object emits no warnings", {
  o <- fit.continuous()
  p <- predict(o$fit, x = o$d$features, tm = o$d$time, id = o$d$id, y = o$d$y)
  on.null.device(expect_no_warning(plot(p)))
})

test_that("print reports the model and returns its argument invisibly", {
  o <- fit.continuous()
  expect_output(print(o$fit), "model")
  expect_output(print(o$fit), "n\\b|sample size|subjects")
  expect_invisible(print(o$fit))
})

test_that("print works on a predict object and on a non-cross-validated fit", {
  o <- fit.continuous()
  p <- predict(o$fit, x = o$d$features, tm = o$d$time, id = o$d$id, y = o$d$y)

  expect_output(print(p), "model")
  expect_output(print(fit.continuous(cv.flag = FALSE)$fit), "model")
})

test_that("display helpers reject objects that are neither grow nor predict", {
  expect_error(print.boostmtree(list(a = 1)))
  expect_error(plot.boostmtree(list(a = 1)))
})

test_that("plot and print cover the non-continuous families", {
  d <- local({
    set.seed(8)
    simLong(
      n = 24, n.time = 3, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  cases <- list(
    binary = function(n) { set.seed(9); rbinom(n, 1, 0.4) },
    nominal = function(n) {
      set.seed(9); factor(sample(c("a", "b", "c"), n, replace = TRUE))
    }
  )
  on.null.device(for (family in names(cases)) {
    # Build the response before seeding the fit. The generators seed
    # themselves, and R forces arguments inside the call, so generating y in
    # the argument list would reset the RNG after set.seed() and leave the fit
    # seed doing nothing.
    y <- cases[[family]](length(d$y))
    set.seed(7)
    fit <- boostmtree(
      d$features, d$time, d$id, y,
      family = family, M = 8, cv.flag = TRUE, verbose = FALSE
    )
    expect_no_warning(plot(fit))
    expect_output(print(fit), "model")
  })
})

test_that("boostmtree.news runs", {
  expect_no_error(suppressWarnings(capture.output(boostmtree.news())))
})
