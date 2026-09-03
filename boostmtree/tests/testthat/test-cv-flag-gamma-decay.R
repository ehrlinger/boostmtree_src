# Regression test for the cv.flag = TRUE in-sample mu freeze.
#
# When the in-sample l.pred/mu refresh is skipped under cv.flag = TRUE, the
# residuals never shrink, phi inflates, the ridge penalty lambda collapses by
# ~3 orders of magnitude, and every boosting step contributes a same-sized
# gamma. The linear predictor then accumulates M unshrunk increments and
# diverges linearly in M.
#
# Both assertions below test that the boosting actually converges, rather than
# pinning stored numerics, so they stay robust across RNG and platform.

test_that("cv.flag = TRUE still yields a converging in-sample boost", {
  set.seed(202609)
  d <- simLong(
    n = 60, n.time = 5, rho = 0.8, model = 2,
    family = "continuous", q = 0
  )$data.list

  fit <- boostmtree(
    d$features, d$time, d$id, d$y,
    family = "continuous", M = 100, cv.flag = TRUE, verbose = FALSE
  )

  # 1. Terminal-node coefficients must decay as the residual shrinks.
  #    Measured: 0.204 when correct, 1.133 when frozen.
  gamma.q <- fit$gamma[[1]]
  node.max <- function(m) max(abs(as.matrix(gamma.q[[m]])[, -1, drop = FALSE]))
  decay.ratio <- node.max(length(gamma.q)) / node.max(1)

  expect_lt(decay.ratio, 0.40)

  # 2. In-sample predictions must stay near the observed response range.
  #    The non-CV path is what predict(), partial.plot() and vimp() read.
  observed <- range(d$y)
  width <- diff(observed)
  lower <- observed[1] - 0.25 * width
  upper <- observed[2] + 0.25 * width

  predicted.mu <- unlist(predict(fit, use.cv.flag = FALSE)$mu)
  expect_true(length(predicted.mu) > 0)
  expect_true(all(is.finite(predicted.mu)))

  predicted <- range(predicted.mu)

  expect_gte(predicted[1], lower)
  expect_lte(predicted[2], upper)
})
