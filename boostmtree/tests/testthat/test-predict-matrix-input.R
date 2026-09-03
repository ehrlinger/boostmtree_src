# Characterization test for matrix versus data.frame covariate input.
#
# `boostmtree()` accepts a matrix for `x`, so a fit trained on a matrix should
# round-trip through `predict(..., x =, tm =, id =, y =)`. An earlier layout
# failed here for nominal fits: the internal predict call rejected the test
# covariates as not matching the training ones, while the data.frame-trained
# equivalent worked. The asymmetry was silent and undocumented.

make.multiclass <- function() {
  set.seed(42)
  n.subject <- 12
  n.obs <- 3
  n <- n.subject * n.obs
  list(
    id = rep(seq_len(n.subject), each = n.obs),
    tm = rep(seq_len(n.obs) - 1, times = n.subject),
    x = matrix(
      rnorm(n * 4), nrow = n,
      dimnames = list(NULL, paste0("x", 1:4))
    ),
    y = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
}

fit.nominal <- function(x, d) {
  set.seed(5)
  boostmtree(
    x = x, tm = d$tm, id = d$id, y = d$y,
    family = "nominal", M = 8, verbose = FALSE
  )
}

test_that("a matrix-trained nominal fit round-trips through predict with newdata", {
  d <- make.multiclass()
  fit <- fit.nominal(d$x, d)

  expect_no_error(predict(fit))
  expect_no_error(predict(fit, x = d$x, tm = d$tm, id = d$id, y = d$y))
})

test_that("matrix and data.frame input give the same nominal fit", {
  d <- make.multiclass()

  from.matrix <- fit.nominal(d$x, d)
  from.frame <- fit.nominal(as.data.frame(d$x), d)

  expect_equal(unlist(from.matrix$mu), unlist(from.frame$mu))
})

test_that("matrix input round-trips for a continuous fit too", {
  d <- make.multiclass()
  set.seed(5)
  y.numeric <- rnorm(length(d$y))
  fit <- boostmtree(
    x = d$x, tm = d$tm, id = d$id, y = y.numeric,
    family = "continuous", M = 8, verbose = FALSE
  )

  expect_no_error(predict(fit, x = d$x, tm = d$tm, id = d$id, y = y.numeric))
})
