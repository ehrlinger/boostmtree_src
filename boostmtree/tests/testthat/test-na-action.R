# `na.action` is documented as accepting "na.omit" or "na.impute", but the
# formal default is written `c("na.omit", "na.impute")[2]`, which evaluates to
# the single string "na.impute". A bare `match.arg()` takes its choices from
# that evaluated default, so "na.omit" was rejected and the documented option
# was unreachable. The `na.omit` branch in
# `boostmtree.remove.missing.subjects()` was therefore dead code.

test_that("both documented na.action values are accepted", {
  d <- local({
    set.seed(4)
    simLong(
      n = 15, n.time = 3, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  fit.with <- function(na) {
    set.seed(11)
    boostmtree(
      d$features, d$time, d$id, d$y,
      family = "continuous", M = 3, verbose = FALSE, na.action = na
    )
  }

  expect_no_error(fit.with("na.impute"))
  expect_no_error(fit.with("na.omit"))
  expect_equal(fit.with("na.omit")$na.action, "na.omit")
  expect_equal(fit.with("na.impute")$na.action, "na.impute")
})

test_that("na.action still defaults to na.impute", {
  d <- local({
    set.seed(4)
    simLong(
      n = 15, n.time = 3, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  set.seed(11)
  fit <- boostmtree(
    d$features, d$time, d$id, d$y,
    family = "continuous", M = 3, verbose = FALSE
  )
  expect_equal(fit$na.action, "na.impute")
})

test_that("an unsupported na.action is still rejected", {
  d <- local({
    set.seed(4)
    simLong(
      n = 15, n.time = 3, rho = 0.8, model = 2,
      family = "continuous", q = 0
    )$data.list
  })
  expect_error(
    boostmtree(
      d$features, d$time, d$id, d$y,
      family = "continuous", M = 3, verbose = FALSE, na.action = "na.fail"
    ),
    "should be one of|'arg' should be"
  )
})
