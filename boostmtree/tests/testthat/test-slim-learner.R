# A fit stores one rfsrc tree per boosting step. Since 2.0.9 each stored tree
# drops what no later call reads (the out-of-bag output, the top-level aliases
# of the data, and the covariate frame that every step shares), and the
# cross-validated coefficients of each left-out subject are one matrix rather
# than a two-element list. These tests pin that the storage changed and the
# results did not. tools/baseline/slim-learner-equivalence.R compares every
# user-facing output bit for bit against a build without the change.

slim.data <- function(n = 30) {
  set.seed(4)
  d <- simLong(n = n, n.time = 3, rho = 0.8, model = 2, family = "continuous", q = 3)$data.list
  x <- d$features
  x$grp <- factor(sample(c("a", "b", "c"), nrow(x), replace = TRUE))
  rows <- split(seq_len(nrow(x)), d$id)
  x$x1[rows[[2]]] <- NA
  x$grp[rows[[5]]] <- NA
  d$features <- x
  d
}

slim.fit <- function(cv.flag = TRUE, na.action = "na.impute") {
  d <- slim.data()
  set.seed(3)
  list(
    fit = boostmtree(d$features, d$time, d$id, d$y, family = "continuous", M = 6,
                     cv.flag = cv.flag, na.action = na.action, mtry = 2, verbose = FALSE),
    d = d
  )
}

# Undo the slimming: what a fit made before 2.0.9 looked like to predict().
unslim <- function(fit) {
  fit$base.learner <- lapply(fit$base.learner, function(by.m) {
    lapply(by.m, boostmtree:::boostmtree.restore.learner, shared.xvar = fit$base.learner.xvar)
  })
  fit$base.learner.xvar <- NULL
  fit
}

test_that("stored trees drop the dead fields and share one covariate frame", {
  fit <- slim.fit()$fit
  expect_s3_class(fit$base.learner.xvar, "data.frame")
  for (learner in fit$base.learner[[1]]) {
    expect_s3_class(learner, "rfsrc")
    expect_false(any(c("xvar", "yvar", "regrOutput", "membership") %in% names(learner)))
    expect_true("xvar" %in% names(learner$forest))
    expect_null(learner$forest$xvar)
    expect_false(is.null(learner$forest$yvar))
    expect_false(is.null(learner$inbag))
  }
})

test_that("a slimmed tree, restored, assigns the same terminal nodes as the full tree", {
  d <- slim.data()
  x <- d$features[!duplicated(d$id), ]
  data <- data.frame(y1 = rnorm(nrow(x)), y2 = rnorm(nrow(x)), x, check.names = FALSE)
  for (na.action in c("na.impute", "na.omit")) {
    tree <- randomForestSRC::rfsrc(Multivar(y1, y2) ~ ., data = data, ntree = 1, mtry = 2,
                                   nodesize = 3, membership = TRUE, na.action = na.action,
                                   nimpute = 1, seed = -7)
    shared <- tree$forest$xvar
    slim <- boostmtree:::boostmtree.slim.learner(tree, shared)
    expect_null(slim$forest$xvar)
    expect_lt(as.numeric(object.size(slim)), as.numeric(object.size(tree)))
    restored <- boostmtree:::boostmtree.restore.learner(slim, shared)
    full <- predict(tree, newdata = x, membership = TRUE, ptn.count = 4,
                    na.action = na.action, importance = "none")$ptn.membership
    again <- predict(restored, newdata = x, membership = TRUE, ptn.count = 4,
                     na.action = na.action, importance = "none")$ptn.membership
    expect_identical(again, full)
  }
})

test_that("a tree whose covariate frame differs keeps its own copy", {
  d <- slim.data()
  x <- d$features[!duplicated(d$id), ]
  data <- data.frame(y1 = rnorm(nrow(x)), x, check.names = FALSE)
  tree <- randomForestSRC::rfsrc(y1 ~ ., data = data, ntree = 1, seed = -7)
  other <- tree$forest$xvar
  other[1, 1] <- other[1, 1] + 1
  slim <- boostmtree:::boostmtree.slim.learner(tree, other)
  expect_identical(slim$forest$xvar, tree$forest$xvar)
  expect_identical(boostmtree:::boostmtree.restore.learner(slim, other)$forest$xvar, tree$forest$xvar)
})

test_that("predict() and vimp() match a fit holding full trees", {
  for (na.action in c("na.impute", "na.omit")) {
    o <- slim.fit(na.action = na.action)
    old <- unslim(o$fit)
    expect_identical(predict(o$fit)$muhat, predict(old)$muhat)
    expect_identical(predict(o$fit, use.cv.flag = TRUE)$muhat, predict(old, use.cv.flag = TRUE)$muhat)
    rows <- o$d$id %in% unique(o$d$id)[1:10]
    p.new <- predict(o$fit, x = o$d$features[rows, ], tm = o$d$time[rows], id = o$d$id[rows], y = o$d$y[rows])
    p.old <- predict(old, x = o$d$features[rows, ], tm = o$d$time[rows], id = o$d$id[rows], y = o$d$y[rows])
    expect_identical(p.new$muhat, p.old$muhat)
    expect_identical(p.new$err.rate, p.old$err.rate)
    set.seed(9)
    v.new <- vimp.boostmtree(o$fit)
    set.seed(9)
    v.old <- vimp.boostmtree(old)
    expect_identical(v.new, v.old)
    set.seed(9)
    v.new <- vimp.boostmtree(p.new)
    set.seed(9)
    v.old <- vimp.boostmtree(p.old)
    expect_identical(v.new, v.old)
  }
})

test_that("cross-validated coefficients are stored one matrix per left-out subject", {
  fit <- slim.fit()$fit
  steps <- Filter(Negate(is.null), unlist(fit$gamma.i.list, recursive = FALSE)[[1]])
  expect_gt(length(steps), 0L)
  for (step in steps) {
    expect_true(is.matrix(step))
    expect_type(step, "double")
  }
})

test_that("compacting a coefficient step keeps the values and label matching", {
  step <- list(node.label = c(2L, 5L, 9L), coefficients = matrix(c(0.1, 0.2, 0.3, 1, 2, 3), 3))
  compact <- boostmtree:::boostmtree.compact.gamma.step(step)
  expect_true(is.matrix(compact))
  expect_lt(as.numeric(object.size(compact)), as.numeric(object.size(step)))
  for (label in c(2L, 5L, 9L)) {
    expect_identical(
      boostmtree:::boostmtree.predict.match.gamma(compact, label),
      boostmtree:::boostmtree.predict.match.gamma(step, label)
    )
    expect_identical(
      boostmtree:::boostmtree.vimp.extract.coefficients(compact, label),
      boostmtree:::boostmtree.vimp.extract.coefficients(step, label)
    )
  }
})

test_that("a stored tree cannot reach randomForestSRC's native predict", {
  fit <- slim.fit(cv.flag = FALSE)$fit
  tree <- fit$base.learner[[1]][[1]]
  expect_identical(class(tree)[1], "boostmtree.learner")
  expect_error(predict(tree), "stored by boostmtree without its training covariates")
  restored <- boostmtree:::boostmtree.restore.learner(tree, fit$base.learner.xvar)
  expect_false(inherits(restored, "boostmtree.learner"))
  expect_false(is.null(restored$forest$xvar))
  # Without boostmtree loaded, predict() dispatches straight to
  # randomForestSRC. Before 2.0.9's class guard this was a segmentation fault,
  # which would take the test runner down, so it runs in a separate R process.
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(tree, path)
  code <- sprintf(paste(
    "suppressMessages(library(randomForestSRC))",
    "tree <- readRDS('%s')",
    "r <- tryCatch(predict(tree), error = function(e) 'error')",
    "cat(if (identical(r, 'error')) 'ERROR' else 'PREDICTED')",
    sep = "; "
  ), path)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(code)),
    stdout = TRUE, stderr = FALSE
  ))
  expect_null(attr(out, "status"))
  expect_identical(tail(out, 1), "ERROR")
})

test_that("compacted node labels match the same way under any scipen", {
  step <- list(node.label = c(1L, 7L, 100000L), coefficients = matrix(1:6 / 10, 3))
  saved <- options(scipen = 999)
  on.exit(options(saved), add = TRUE)
  compact <- boostmtree:::boostmtree.compact.gamma.step(step)
  expect_true(is.matrix(compact))
  for (scipen in c(0, 999, -5)) {
    options(scipen = scipen)
    for (label in c(1L, 7L, 100000L)) {
      expect_identical(
        boostmtree:::boostmtree.predict.match.gamma(compact, label),
        boostmtree:::boostmtree.predict.match.gamma(step, label)
      )
      expect_identical(
        boostmtree:::boostmtree.vimp.extract.coefficients(compact, label),
        boostmtree:::boostmtree.vimp.extract.coefficients(step, label)
      )
    }
  }
  fractional <- list(node.label = c(1.5, 2), coefficients = matrix(1:4 / 10, 2))
  expect_identical(boostmtree:::boostmtree.compact.gamma.step(fractional), fractional)
})
