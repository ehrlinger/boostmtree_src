# progress = reports how far a fit or a partial.plot() has got, for callers that
# cannot see verbose's progress bar (knitr captures it into the report).

progress.data <- function() {
  set.seed(2)
  simLong(n = 30, n.time = 4, rho = 0.8, model = 2, family = "continuous", q = 0)$data.list
}

test_that("progress is called once per boosting iteration, in order, with and without cv.flag", {
  d <- progress.data()
  for (cv in c(FALSE, TRUE)) {
    seen <- list()
    set.seed(7)
    boostmtree(d$features, d$time, d$id, d$y, family = "continuous", M = 10, cv.flag = cv, verbose = FALSE,
               progress = function(step, steps) seen[[length(seen) + 1L]] <<- c(step, steps))
    expect_equal(seen, lapply(1:10, function(m) c(m, 10)))
  }
})

test_that("a failing progress callback warns once and leaves the fit unchanged", {
  d <- progress.data()
  set.seed(7)
  a <- boostmtree(d$features, d$time, d$id, d$y, family = "continuous", M = 10, verbose = FALSE)
  calls <- 0L
  set.seed(7)
  expect_warning(
    b <- boostmtree(d$features, d$time, d$id, d$y, family = "continuous", M = 10, verbose = FALSE,
                    progress = function(step, steps) {
                      calls <<- calls + 1L
                      stop("boom")
                    }),
    "boom"
  )
  expect_identical(calls, 1L)
  expect_setequal(names(b), names(a))
  # the forests in base.learner carry wall-clock timings, so they cannot be compared
  keep <- setdiff(names(a), c("call", "base.learner"))
  expect_equal(unclass(b)[keep], unclass(a)[keep])
})

test_that("progress must be NULL or a function", {
  d <- progress.data()
  expect_error(boostmtree(d$features, d$time, d$id, d$y, family = "continuous", M = 2, verbose = FALSE, progress = 1),
               "`progress` must be NULL or a function")
})

test_that("partial.plot reports one step per (variable, grid point), with steps fixed from the start", {
  d <- progress.data()
  set.seed(7)
  fit <- boostmtree(d$features, d$time, d$id, d$y, family = "continuous", M = 5, cv.flag = TRUE, verbose = FALSE)
  seen <- list()
  partial.plot(fit, x.var.names = names(d$features)[1:2], n.points = 3, output = "data", verbose = FALSE,
               progress = function(step, steps) seen[[length(seen) + 1L]] <<- c(step, steps))
  steps <- vapply(seen, `[[`, 0, 2)
  expect_true(length(seen) >= 2L && length(seen) <= 6L)
  expect_true(all(steps == length(seen)))
  expect_equal(vapply(seen, `[[`, 0, 1), seq_along(seen))
})

test_that("partial.plot's result does not depend on progress", {
  d <- progress.data()
  set.seed(7)
  fit <- boostmtree(d$features, d$time, d$id, d$y, family = "continuous", M = 5, verbose = FALSE)
  a <- partial.plot(fit, n.points = 3, output = "data", verbose = FALSE)
  b <- partial.plot(fit, n.points = 3, output = "data", verbose = FALSE, progress = function(step, steps) NULL)
  expect_equal(b[setdiff(names(b), "call")], a[setdiff(names(a), "call")])
})
