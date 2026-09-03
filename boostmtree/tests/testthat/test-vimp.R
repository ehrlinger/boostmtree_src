# Coverage for vimp.boostmtree(), the largest uncovered file in the package
# (362 lines at 0%). The multinomial accumulation path, the predict-object
# path that the documentation examples use, and the joint variant were all
# entirely unexercised.

fit.continuous <- function(cv.flag = TRUE) {
  d <- local({
    set.seed(2)
    simLong(
      n = 26, n.time = 4, rho = 0.8, model = 2,
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

fit.categorical <- function(family, y.fun, cv.flag = TRUE) {
  d <- local({
    set.seed(8)
    simLong(
      n = 24, n.time = 3, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  set.seed(7)
  list(
    fit = boostmtree(
      d$features, d$time, d$id, y.fun(length(d$y)),
      family = family, M = 8, cv.flag = cv.flag, verbose = FALSE
    ),
    d = d
  )
}

expect.well.formed.vimp <- function(v, n.x) {
  expect_s3_class(v, "vimp.boostmtree")
  expect_true(is.matrix(v$main))
  expect_true(is.matrix(v$interaction))
  expect_equal(nrow(v$main), n.x)
  expect_equal(dim(v$main), dim(v$interaction))
  expect_true(all(is.finite(v$main)))
  expect_true(all(is.finite(v$interaction)))
  expect_true(all(is.finite(v$time.effect)))
  # multi-class families carry one baseline per class, so this is a vector
  expect_true(all(is.finite(v$baseline.rmse)))
}

test_that("vimp on a grow object is well formed for a continuous response", {
  o <- fit.continuous()
  v <- vimp.boostmtree(o$fit)

  expect.well.formed.vimp(v, length(o$fit$x.var.names))
  expect_equal(rownames(v$main), o$fit$x.var.names)
  expect_equal(v$source, "grow")
  expect_false(v$joint)
})

test_that("vimp works on a predict object, the form the documentation uses", {
  o <- fit.continuous()
  p <- predict(o$fit, x = o$d$features, tm = o$d$time, id = o$d$id, y = o$d$y)
  v <- vimp.boostmtree(p)

  expect.well.formed.vimp(v, length(o$fit$x.var.names))
  expect_equal(v$source, "predict")
})

test_that("x.names restricts vimp to the named covariates", {
  o <- fit.continuous()
  v <- vimp.boostmtree(o$fit, x.names = c("x1", "x3"))

  expect_equal(rownames(v$main), c("x1", "x3"))
  expect_equal(nrow(v$main), 2L)
})

test_that("joint vimp returns a single combined row rather than one per covariate", {
  # Regression: with the default x.names = NULL, the joint collapse was skipped
  # while the row count still became 1, so labelling the matrix failed with
  # "length of 'dimnames' [1] not equal to array extent". Supplying x.names
  # took the other branch and masked it.
  o <- fit.continuous()
  single <- vimp.boostmtree(o$fit, joint = FALSE)
  joint <- vimp.boostmtree(o$fit, joint = TRUE)

  expect_true(joint$joint)
  expect_equal(nrow(joint$main), 1L)
  expect_equal(rownames(joint$main), "joint.vimp")
  expect_lt(nrow(joint$main), nrow(single$main))
  expect_true(all(is.finite(joint$main)))
})

test_that("joint vimp also works on a predict object and with x.names supplied", {
  o <- fit.continuous()
  p <- predict(o$fit, x = o$d$features, tm = o$d$time, id = o$d$id, y = o$d$y)

  expect_equal(rownames(vimp.boostmtree(p, joint = TRUE)$main), "joint.vimp")
  expect_equal(
    rownames(vimp.boostmtree(o$fit, x.names = c("x1", "x2"), joint = TRUE)$main),
    "joint.vimp"
  )
})

test_that("vimp covers the binary family", {
  o <- fit.categorical("binary", function(n) {
    set.seed(9); rbinom(n, 1, 0.4)
  })
  v <- vimp.boostmtree(o$fit)
  expect.well.formed.vimp(v, length(o$fit$x.var.names))
  expect_equal(v$family, "binary")
})

test_that("vimp covers the nominal family and accumulates over classes", {
  o <- fit.categorical("nominal", function(n) {
    set.seed(9); factor(sample(c("a", "b", "c"), n, replace = TRUE))
  })
  v <- vimp.boostmtree(o$fit)

  expect.well.formed.vimp(v, length(o$fit$x.var.names))
  expect_equal(v$family, "nominal")
  expect_gt(ncol(v$main), 1L)
})

test_that("vimp covers the ordinal family", {
  o <- fit.categorical("ordinal", function(n) {
    set.seed(9)
    factor(sample(c("low", "mid", "high"), n, replace = TRUE),
           levels = c("low", "mid", "high"), ordered = TRUE)
  })
  v <- vimp.boostmtree(o$fit)

  expect.well.formed.vimp(v, length(o$fit$x.var.names))
  expect_equal(v$family, "ordinal")
})

test_that("grow-object vimp requires a cross-validated fit and says so", {
  # Documented constraint: the grow path measures OOB RMSE, which only exists
  # when the fit carried cv.flag = TRUE.
  o <- fit.continuous(cv.flag = FALSE)
  expect_error(vimp.boostmtree(o$fit), "requires `cv.flag = TRUE`")
})

test_that("a non-cross-validated fit can still be scored through a predict object", {
  o <- fit.continuous(cv.flag = FALSE)
  p <- predict(o$fit, x = o$d$features, tm = o$d$time, id = o$d$id, y = o$d$y)
  expect.well.formed.vimp(vimp.boostmtree(p), length(o$fit$x.var.names))
})

test_that("vimp rejects objects it cannot handle", {
  expect_error(vimp.boostmtree(list(a = 1)), "only works for")
  expect_error(vimp.boostmtree(42), "only works for")
})

test_that("plotting a vimp object emits no warnings and returns invisibly", {
  o <- fit.continuous()
  v <- vimp.boostmtree(o$fit)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_warning(plot(v))
  expect_no_warning(plot(vimp.boostmtree(o$fit, joint = TRUE)))
})
