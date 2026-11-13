# Characterization test for the NA-removal path in predict().
#
# An earlier layout subset the covariate frame on that path with a missing
# comma (`X[index, drop = FALSE]`), which selects columns rather than rows.
# The neighbouring `tm`, `id` and `y` were subset correctly by the same index,
# so the covariates silently fell out of step with them and predictions were
# corrupted rather than erroring.
#
# The invariant below is what "in step" means, and it is what a wrong-dimension
# subset breaks: dropping a subject by making its covariates NA must give the
# surviving subjects exactly the same predictions as never supplying that
# subject at all.

make.longitudinal <- function(n = 25, n.time = 4) {
  set.seed(4)
  simLong(
    n = n, n.time = n.time, rho = 0.8, model = 2,
    family = "continuous", q = 0
  )$data.list
}

fit.omit <- function(d) {
  set.seed(11)
  boostmtree(
    d$features, d$time, d$id, d$y,
    family = "continuous", M = 10, cv.flag = FALSE,
    verbose = FALSE, na.action = "na.omit"
  )
}

test_that("a subject dropped via NA leaves the other subjects' predictions unchanged", {
  d <- make.longitudinal()
  fit <- fit.omit(d)

  dropped.id <- unique(d$id)[1L]
  is.dropped <- d$id == dropped.id

  # (a) same data, but the dropped subject's covariates are NA
  x.na <- d$features
  x.na[is.dropped, 1L] <- NA_real_
  with.na <- predict(fit, x = x.na, tm = d$time, id = d$id, y = d$y)

  # (b) same data, but the dropped subject's rows are absent entirely
  without <- predict(
    fit,
    x = d$features[!is.dropped, , drop = FALSE],
    tm = d$time[!is.dropped],
    id = d$id[!is.dropped],
    y = d$y[!is.dropped]
  )

  expect_false(dropped.id %in% with.na$id.unique)
  expect_equal(with.na$id.unique, without$id.unique)
  expect_equal(with.na$mu, without$mu)
})

test_that("the NA-removal path keeps covariates aligned with tm, id and y", {
  d <- make.longitudinal()
  fit <- fit.omit(d)

  dropped.id <- unique(d$id)[2L]
  is.dropped <- d$id == dropped.id

  x.na <- d$features
  x.na[is.dropped, 1L] <- NA_real_
  p <- predict(fit, x = x.na, tm = d$time, id = d$id, y = d$y)

  expect_equal(nrow(p$x), length(p$id.unique))
  expect_equal(ncol(p$x), ncol(d$features))
  expect_equal(colnames(p$x), colnames(d$features))
  expect_equal(length(p$mu), length(p$id.unique))
  expect_equal(vapply(p$time, length, integer(1)), p$ni)
  expect_true(all(is.finite(unlist(p$mu))))
})

# With `x` alone (no `tm` and `id`), and with `partial = TRUE`, an NA row used
# to reach randomForestSRC, which dropped it, so the terminal-node lookup no
# longer lined up with the rows and predict() stopped with "subscript out of
# bounds". The row is now dropped first, with a warning, and `id.unique` names
# the rows of `x` that were kept.
test_that("an na.omit fit drops NA rows of x given without tm and id", {
  d <- make.longitudinal()
  fit <- fit.omit(d)
  x <- d$features[!duplicated(d$id), , drop = FALSE][1:6, , drop = FALSE]
  x.na <- x
  x.na[2L, 1L] <- NA_real_
  clean <- predict(fit, x = x[-2L, , drop = FALSE])
  expect_warning(p <- predict(fit, x = x.na), "1 row\\(s\\) of `x`")
  expect_identical(p$id.unique, c(1L, 3:6))
  expect_identical(p$muhat, clean$muhat)
  tm <- sort(unique(d$time))
  clean <- predict(fit, x = x[-2L, , drop = FALSE], tm = tm, partial = TRUE)
  expect_warning(p <- predict(fit, x = x.na, tm = tm, partial = TRUE), "1 row\\(s\\) of `x`")
  expect_identical(p$id.unique, c(1L, 3:6))
  expect_identical(p$muhat, clean$muhat)
  x.na[, 1L] <- NA_real_
  expect_error(predict(fit, x = x.na), "Every row of `x` has a missing covariate")
})

# The cross-validated coefficients are looked up by fitted subject position.
# One extra row plus one dropped NA row leaves the count equal to the fitted
# subject count, but the rows are no longer the fitted set in order, so
# use.cv.flag must fall back to the ordinary coefficients.
test_that("a dropped NA row turns off use.cv.flag for partial predictions", {
  d <- make.longitudinal(n = 20, n.time = 3)
  set.seed(11)
  fit <- boostmtree(d$features, d$time, d$id, d$y, family = "continuous", M = 5,
                    cv.flag = TRUE, verbose = FALSE, na.action = "na.omit")
  tm <- sort(unique(d$time))
  x <- rbind(fit$x, fit$x[1L, , drop = FALSE])
  x[2L, 1L] <- NA_real_
  warnings <- character(0)
  p <- withCallingHandlers(
    predict(fit, x = x, tm = tm, partial = TRUE, use.cv.flag = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("using the ordinary fitted coefficients", warnings)))
  expect_identical(p$muhat, suppressWarnings(predict(fit, x = x, tm = tm, partial = TRUE))$muhat)
})
