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
  if (.Platform$OS.type == "windows") {
    expect_warning(
      ctl <- boostmtree:::boostmtree.normalize.control(boostmtree.control(cv.cores = 3)),
      "Windows"
    )
    expect_identical(ctl$cv.cores, 1L)
  } else {
    expect_identical(boostmtree:::boostmtree.normalize.control(boostmtree.control(cv.cores = 3))$cv.cores, 3L)
  }
})

# mclapply() returns NULL for a worker that dies without delivering, rather
# than a try-error. Assigning NULL into the state list deletes the subject and
# the fit carries on one subject short, so a NULL must stop the fit.
test_that("an undelivered worker result stops the fit", {
  skip_on_os("windows")
  local_mocked_bindings(
    mclapply = function(X, FUN, ...) {
      out <- lapply(X, FUN)
      out[1] <- list(NULL)
      out
    },
    .package = "parallel"
  )
  set.seed(2026)
  d <- simLong(n = 20, n.time = 4, model = 2, family = "continuous")$data.list
  expect_error(
    boostmtree(
      d$features, d$time, d$id, d$y,
      family = "continuous", M = 3, cv.flag = TRUE, verbose = FALSE,
      control = boostmtree.control(cv.cores = 2L)
    ),
    "not delivered"
  )
})

# The ordinal cross-validated update raises each threshold's state to the one
# below it, element by element. Neither test above reaches that branch: it runs
# only for an ordinal response with cv.flag = TRUE. The state is internal, and
# the cross-validated mu is rescaled per threshold, so mu itself need not be
# ordered; the branch is therefore checked on the step directly.
test_that("the ordinal CV step raises each subject's state to the threshold below", {
  ni <- c(2L, 3L, 1L)
  n <- length(ni)
  rows <- unname(split(seq_len(sum(ni)), rep.int(seq_len(n), ni)))
  set.seed(1)
  state <- lapply(seq_len(n), function(i) rnorm(sum(ni)))
  below <- lapply(seq_len(n), function(i) rnorm(sum(ni)))
  design <- lapply(ni, function(k) cbind(1, seq_len(k)))
  step <- boostmtree:::boostmtree.update.cv.step(
    cv.state.q = state, subject.rows = rows,
    rfsrc.obj = list(inbag = rep(1L, n)),
    membership = rep(1, n), membership.original = rep(1, n),
    time.design = design, nu.vec = c(1, 1), lambda = 1,
    y.current = lapply(ni, numeric), inv.v.matrix = lapply(ni, diag),
    family = "ordinal", nr.iter = 1L, pen.lsq.matrix = diag(2),
    previous.ordinal = below
  )
  expect_identical(step$l.pred.db.i.q, Map(pmax, state, below))
})

test_that("an ordinal cross-validated fit is core-invariant", {
  skip_on_cran()
  skip_on_os("windows")
  fit.with <- function(cores) {
    set.seed(2026)
    d <- simLong(n = 40, n.time = 5, model = 2, family = "continuous")$data.list
    y <- cut(d$y, quantile(d$y, c(0, 0.6, 0.85, 1)), labels = 0:2,
             include.lowest = TRUE, ordered_result = TRUE)
    set.seed(77)
    boostmtree(
      d$features, d$time, d$id, y,
      family = "ordinal", M = 8, nu = 0.5, cv.flag = TRUE, verbose = FALSE,
      control = boostmtree.control(cv.cores = cores)
    )
  }
  serial <- fit.with(1L)
  forked <- fit.with(2L)
  expect_length(serial$mu, 2L)
  expect_true(all(is.finite(unlist(serial$err.rate))))
  expect_identical(forked$mu, serial$mu)
  expect_identical(forked$err.rate, serial$err.rate)
  expect_identical(forked$m.opt, serial$m.opt)
})

# A multithreaded BLAS starts one thread per CPU in every forked worker, so the
# parallel out-of-bag update must hold each worker, and the parent while it
# waits, to one BLAS thread. The BLAS here is mocked: each mocked worker starts
# at 8 threads, as a fresh fork can, and mclapply() runs in this process so the
# calls are visible.
test_that("the parallel CV step holds each worker and the parent to one BLAS thread", {
  skip_on_os("windows")
  skip_on_cran()
  blas <- new.env()
  blas$threads <- 8L
  blas$set <- integer(0)
  blas$worker.calls <- 0L
  local_mocked_bindings(
    blas_get_num_procs = function() blas$threads,
    blas_set_num_threads = function(threads) {
      blas$set <- c(blas$set, threads)
      blas$threads <- as.integer(threads)
      invisible(NULL)
    },
    .package = "RhpcBLASctl"
  )
  local_mocked_bindings(
    mclapply = function(X, FUN, ...) {
      lapply(X, function(x) {
        blas$threads <- 8L
        blas$worker.calls <- blas$worker.calls + 1L
        FUN(x)
      })
    },
    .package = "parallel"
  )
  set.seed(2026)
  d <- simLong(n = 20, n.time = 4, model = 2, family = "continuous")$data.list
  fit <- function(cores) {
    boostmtree(
      d$features, d$time, d$id, d$y,
      family = "continuous", M = 2, cv.flag = TRUE, verbose = FALSE,
      control = boostmtree.control(cv.cores = cores)
    )
  }
  fit(2L)
  steps <- 2L
  expect_gt(blas$worker.calls, 0L)
  # One parent drop and one restore per step, and one drop per worker call.
  expect_identical(sum(blas$set == 1L), blas$worker.calls + steps)
  expect_identical(sum(blas$set == 8L), steps)
  expect_identical(blas$threads, 8L)
  blas$set <- integer(0)
  fit(1L)
  expect_length(blas$set, 0L)
})
