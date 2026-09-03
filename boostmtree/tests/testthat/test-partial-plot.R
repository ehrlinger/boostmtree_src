# Coverage for partial.plot(), which had none despite being a main entry point
# and the place a real analysis bug surfaced: a subsetted partial plot silently
# switched coefficient paths and produced values far outside the response range.

make.fit <- function(cv.flag = TRUE) {
  d <- local({
    set.seed(2)
    simLong(
      n = 30, n.time = 4, rho = 0.8, model = 2,
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

test_that("output = 'data' returns one curve frame per covariate, shaped by n.points and time.points", {
  o <- make.fit()
  p <- partial.plot(
    o$fit, x.var.names = c("x1", "x2"),
    time.points = c(0.5, 2), n.points = 5, output = "data"
  )

  expect_s3_class(p, "partial.plot.boostmtree")
  expect_named(p$curves, c("x1", "x2"))
  for (nm in c("x1", "x2")) {
    expect_s3_class(p$curves[[nm]], "data.frame")
    expect_lte(nrow(p$curves[[nm]]), 5L)
    expect_equal(ncol(p$curves[[nm]]), 1L + length(p$time.points))
    expect_equal(names(p$curves[[nm]])[1], "x")
    expect_true(all(is.finite(as.matrix(p$curves[[nm]][, -1, drop = FALSE]))))
  }
})

test_that("requested time points snap to observed times", {
  o <- make.fit()
  p <- partial.plot(
    o$fit, x.var.names = "x1",
    time.points = c(0.5, 2), n.points = 4, output = "data"
  )
  expect_length(p$time.points, 2L)
  expect_true(all(p$time.points %in% o$fit$time.unique))
})

test_that("subset restricts the averaged subjects and forces the non-CV path with a warning", {
  # This is the behaviour that produced impossible values in a real analysis:
  # supplying `subset` silently downgrades use.cv.flag, because the CV
  # coefficients are subject-indexed and do not survive subsetting.
  o <- make.fit(cv.flag = TRUE)
  keep <- o$fit$x$x1 > stats::median(o$fit$x$x1)

  expect_warning(
    p <- partial.plot(
      o$fit, x.var.names = "x1", subset = keep,
      n.points = 4, output = "data", use.cv.flag = TRUE
    ),
    "only supported when all fitted subjects are used"
  )
  expect_false(p$use.cv.flag)

  # the x grid is drawn from the subset, not the full cohort
  expect_gte(min(p$curves$x1$x), min(o$fit$x$x1[keep]))
})

test_that("subset without use.cv.flag = TRUE needs no warning", {
  o <- make.fit(cv.flag = TRUE)
  keep <- o$fit$x$x1 > stats::median(o$fit$x$x1)
  expect_no_warning(
    partial.plot(
      o$fit, x.var.names = "x1", subset = keep,
      n.points = 4, output = "data", use.cv.flag = FALSE
    )
  )
})

test_that("conditioning on a covariate value succeeds and changes the curve", {
  o <- make.fit()
  base <- partial.plot(
    o$fit, x.var.names = "x1", n.points = 4, output = "data",
    time.points = 1
  )
  cond <- partial.plot(
    o$fit, x.var.names = "x1", n.points = 4, output = "data",
    time.points = 1,
    conditional.x.var.names = "x2",
    conditional.values = max(o$fit$x$x2)
  )

  expect_equal(dim(base$curves$x1), dim(cond$curves$x1))
  expect_false(isTRUE(all.equal(base$curves$x1[[2]], cond$curves$x1[[2]])))
})

test_that("conditioning arguments are validated", {
  o <- make.fit()
  expect_error(
    partial.plot(o$fit, x.var.names = "x1", output = "data",
                 conditional.x.var.names = "x2"),
    "must be supplied together"
  )
  expect_error(
    partial.plot(o$fit, x.var.names = "x1", output = "data",
                 conditional.x.var.names = c("x2", "x3"),
                 conditional.values = 1),
    "same length"
  )
  expect_error(
    partial.plot(o$fit, x.var.names = "x1", output = "data",
                 conditional.x.var.names = "nope", conditional.values = 1),
    "not found"
  )
})

test_that("an unmatched covariate name is rejected", {
  o <- make.fit()
  expect_error(
    partial.plot(o$fit, x.var.names = "not_a_variable", output = "data"),
    "does not match any fitted covariate"
  )
})

test_that("plotting emits no warnings or conditions", {
  # The two v2.0.1 label bugs both slipped through because nothing asserted the
  # absence of warnings around a plotting call.
  o <- make.fit()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(partial.plot(o$fit, x.var.names = "x1", n.points = 4))
})
