# The out-of-bag update forks one refit per left-out subject when
# control$cv.cores > 1. The refits draw no random numbers, so the fit must be
# identical, not merely close, to the sequential one. Compare with identical()
# rather than a tolerance: a tolerance would hide an RNG stream disturbed
# inside the forked loop, which is the failure this guards against.

test_that("cv.cores > 1 gives the identical cross-validated fit", {
  skip_on_os("windows")
  skip_on_cran()
  fit.with <- function(cores) {
    set.seed(2026)
    d <- simLong(n = 30, n.time = 4, model = 2, family = "continuous")$data.list
    set.seed(77)
    fit <- boostmtree(
      d$features, d$time, d$id, d$y,
      family = "continuous", M = 15, nu = 0.5, cv.flag = TRUE, verbose = FALSE,
      control = boostmtree.control(cv.cores = cores)
    )
    list(mu = fit$mu, err.rate = fit$err.rate, m.opt = fit$m.opt, rng = .Random.seed)
  }
  serial <- fit.with(1L)
  forked <- fit.with(2L)
  expect_identical(forked$mu, serial$mu)
  expect_identical(forked$err.rate, serial$err.rate)
  expect_identical(forked$m.opt, serial$m.opt)
  expect_identical(forked$rng, serial$rng)
})

test_that("cv.cores is validated", {
  expect_error(boostmtree:::boostmtree.normalize.control(boostmtree.control(cv.cores = 0)), "cv.cores")
  expect_error(boostmtree:::boostmtree.normalize.control(boostmtree.control(cv.cores = NA)), "cv.cores")
  expect_identical(boostmtree:::boostmtree.normalize.control(boostmtree.control(cv.cores = 3))$cv.cores, 3L)
})
