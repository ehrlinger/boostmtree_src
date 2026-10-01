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

# Categorical families. Until now only the continuous path was exercised, so
# the per-class nesting of `curves` and the probability scale were unchecked.

make.categorical <- function(family, y.fun) {
  d <- local({
    set.seed(8)
    simLong(
      n = 24, n.time = 3, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  # Build the response before seeding the fit (see test-marginal-plot.R).
  y <- y.fun(length(d$y))
  set.seed(7)
  boostmtree(
    d$features, d$time, d$id, y,
    family = family, M = 8, cv.flag = FALSE, verbose = FALSE
  )
}

curve.values <- function(curve) as.matrix(curve[, -1, drop = FALSE])

test_that("partial.plot returns one probability-scale curve set for a binary response", {
  fit <- make.categorical("binary", function(n) {
    set.seed(9); rbinom(n, 1, 0.4)
  })
  expect_no_warning(
    p <- partial.plot(fit, x.var.names = "x1", time.points = 1,
                      n.points = 4, output = "data", verbose = FALSE)
  )

  expect_equal(p$family, "binary")
  expect_length(p$response.labels, 1L)
  expect_named(p$curves, "x1")
  v <- curve.values(p$curves$x1)
  expect_true(all(is.finite(v)))
  expect_true(all(v >= 0 & v <= 1))
})

test_that("partial.plot nests curves by class for a nominal response", {
  fit <- make.categorical("nominal", function(n) {
    set.seed(9); factor(sample(c("a", "b", "c"), n, replace = TRUE))
  })
  expect_no_warning(
    p <- partial.plot(fit, x.var.names = "x1", time.points = c(1, 2),
                      n.points = 4, output = "data", verbose = FALSE)
  )

  expect_gt(length(p$response.labels), 1L)
  expect_named(p$curves, p$response.labels)
  for (per.class in p$curves) {
    v <- curve.values(per.class$x1)
    expect_true(all(is.finite(v)))
    expect_true(all(v >= 0 & v <= 1))
  }
  # The curves are the non-reference classes, so jointly they can take at
  # most the whole probability mass; bounding each one alone does not imply it.
  totals <- Reduce(`+`, lapply(p$curves, function(cl) curve.values(cl$x1)))
  expect_true(all(totals <= 1 + 1e-8))
})

test_that("partial.plot with prob.class = TRUE gives ordinal class probabilities that sum to one", {
  fit <- make.categorical("ordinal", function(n) {
    set.seed(9); factor(sample(1:3, n, replace = TRUE), ordered = TRUE)
  })

  cumulative <- partial.plot(fit, x.var.names = "x1", time.points = c(1, 2),
                             n.points = 4, output = "data", verbose = FALSE)
  expect_length(cumulative$curves, length(cumulative$response.labels))
  levels.v <- lapply(cumulative$curves, function(cl) curve.values(cl$x1))
  for (v in levels.v) {
    expect_true(all(v >= 0 & v <= 1))
  }
  # Cumulative probabilities P(Y <= k) must not cross or reverse across
  # thresholds, at every grid point and time.
  for (k in seq_along(levels.v)[-1]) {
    expect_true(all(levels.v[[k]] >= levels.v[[k - 1]] - 1e-8))
  }

  expect_no_warning(
    by.class <- partial.plot(fit, x.var.names = "x1", time.points = 1,
                             n.points = 4, prob.class = TRUE,
                             output = "data", verbose = FALSE)
  )
  expect_equal(by.class$response.labels, as.character(fit$y.levels))
  expect_length(by.class$curves, length(fit$y.levels))
  totals <- Reduce(`+`, lapply(by.class$curves, function(cl) curve.values(cl$x1)))
  expect_equal(unname(totals), matrix(1, nrow(totals), ncol(totals)),
               tolerance = 1e-8)
})

test_that("prob.class = TRUE is downgraded with a warning outside the ordinal family", {
  fit <- make.categorical("binary", function(n) {
    set.seed(9); rbinom(n, 1, 0.4)
  })
  default <- partial.plot(fit, x.var.names = "x1", time.points = 1,
                          n.points = 4, output = "data", verbose = FALSE)
  expect_warning(
    downgraded <- partial.plot(fit, x.var.names = "x1", time.points = 1,
                               n.points = 4, prob.class = TRUE,
                               output = "data", verbose = FALSE),
    "only used for the ordinal family"
  )
  # Ignored means ignored: the result matches the default call.
  expect_false(downgraded$prob.class)
  expect_identical(downgraded$curves, default$curves)
  expect_identical(downgraded$response.labels, default$response.labels)
})

test_that("partial.plot draws multi-panel plots for nominal and ordinal responses without warnings", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  nominal <- make.categorical("nominal", function(n) {
    set.seed(9); factor(sample(c("a", "b", "c"), n, replace = TRUE))
  })
  ordinal <- make.categorical("ordinal", function(n) {
    set.seed(9); factor(sample(1:3, n, replace = TRUE), ordered = TRUE)
  })

  cases <- list(
    list(fit = nominal, prob.class = FALSE),
    list(fit = ordinal, prob.class = FALSE),
    list(fit = ordinal, prob.class = TRUE)
  )
  for (case in cases) {
    expect_no_warning(
      p <- partial.plot(case$fit, x.var.names = c("x1", "x2"),
                        time.points = 1, n.points = 4,
                        prob.class = case$prob.class, verbose = FALSE)
    )
    expect_s3_class(p, "partial.plot.boostmtree")
    expect_gt(length(p$response.labels), 1L)
    expect_named(p$curves, p$response.labels)
  }
})
