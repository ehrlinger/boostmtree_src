# Regression tests for the step-at-a-time boosting path in predict().
#
# predict(..., partial = TRUE) used to hold every boosting step's linear
# predictor for every subject on the full time grid at once (n.q x M x n
# vectors), twice over: in the main path and in boostmtree.predict.extrapolate().
# On a real ordinal study that was ~55 GB per call, and partial.plot() makes one
# call per grid value. Both paths now keep a running sum. The reference
# implementations below are the old full-path code, kept here so the rolling
# version is checked against it bit for bit.

old.extrapolate <- function(object, membership, gamma, m.opt) {
  n <- length(membership[[1]][[1]])
  n.q <- object$n.q
  nu <- object$nu
  if (length(nu) == 1L) nu <- rep(nu, 2L)
  nu.vec <- c(nu[1], rep(nu[2], ncol(object$x.tm) - 1L))
  x.tm.grid <- object$x.tm
  temp <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n), function(i) {
      lapply(seq_len(m.opt[q]), function(m) {
        gamma.i <- boostmtree:::boostmtree.predict.match.gamma(gamma[[q]][[m]], membership[[q]][[m]][i])
        as.vector(x.tm.grid %*% c(gamma.i[1, ] * nu.vec))
      })
    })
  })
  db <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n), function(i) {
      path <- temp[[q]][[i]]
      if (identical(object$family, "ordinal") && q > 1L) {
        for (m in seq_len(m.opt[q])) path[[m]] <- pmax(path[[m]], temp[[q - 1L]][[i]][[m]])
      }
      Reduce("+", path)
    })
  })
  ref <- if (identical(object$family, "nominal")) {
    lapply(seq_len(n), function(i) {
      log((1 + Reduce("+", lapply(seq_len(n.q), function(q) exp(db[[q]][[i]]))))^(-1))
    })
  } else {
    lapply(seq_len(n), function(i) rep(0, nrow(x.tm.grid)))
  }
  lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n), function(i) {
      boostmtree:::boostmtree.get.mu((db[[q]][[i]] + ref[[i]]) * object$y.sd + object$y.mean, family = object$family)
    })
  })
}

old.mu <- function(object, pred, m.opt) {
  n <- pred$n
  ni <- pred$ni
  n.q <- object$n.q
  m.max <- max(m.opt)
  nu.vec <- pred$nu.vec
  td <- pred$time.design
  membership <- pred$membership
  gamma <- pred$gamma
  db <- lapply(seq_len(n.q), function(q) vector("list", m.max))
  for (q in seq_len(n.q)) {
    for (m in seq_len(m.max)) {
      previous <- if (m == 1L) lapply(seq_len(n), function(i) rep(0, ni[i])) else db[[q]][[m - 1L]]
      if (!pred$use.cv.flag) {
        g <- boostmtree:::boostmtree.predict.match.gamma(gamma[[q]][[m]], membership[[q]][[m]])
      }
      oob <- if (pred$use.cv.flag) which(object$base.learner[[q]][[m]]$inbag == 0) else NULL
      db[[q]][[m]] <- lapply(seq_len(n), function(i) {
        if (!pred$use.cv.flag) {
          updated <- previous[[i]] + c(td[[i]] %*% c(g[i, ] * nu.vec))
        } else {
          updated <- previous[[i]]
          if (i %in% oob) {
            g.i <- boostmtree:::boostmtree.predict.match.gamma(gamma[[q]][[m]][[i]], membership[[q]][[m]][match(i, oob)])
            updated <- updated + c(td[[i]] %*% c(g.i[1, ] * nu.vec))
          }
        }
        if (identical(object$family, "ordinal") && q > 1L) pmax(updated, db[[q - 1L]][[m]][[i]]) else updated
      })
    }
  }
  lapply(seq_len(n.q), function(q) {
    m <- m.opt[q]
    ref <- if (identical(object$family, "nominal")) {
      lapply(seq_len(n), function(i) {
        log((1 + Reduce("+", lapply(seq_len(n.q), function(qq) exp(db[[qq]][[m]][[i]]))))^(-1))
      })
    } else {
      lapply(seq_len(n), function(i) rep(0, ni[i]))
    }
    lapply(seq_len(n), function(i) {
      boostmtree:::boostmtree.get.mu((db[[q]][[m]][[i]] + ref[[i]]) * object$y.sd + object$y.mean, family = object$family)
    })
  })
}

as.q.list <- function(x, family) {
  if (family %in% c("nominal", "ordinal")) x else list(x)
}

make.family.fit <- function(family, cv.flag = FALSE) {
  set.seed(3)
  d <- simLong(
    n = 25, n.time = 4, rho = 0.8, model = 2,
    family = if (identical(family, "binary")) "binary" else "continuous", q = 0
  )$data.list
  y <- d$y
  if (family %in% c("ordinal", "nominal")) {
    y <- cut(y, stats::quantile(y, 0:4 / 4), include.lowest = TRUE, labels = FALSE)
  }
  set.seed(7)
  fit <- boostmtree(d$features, d$time, d$id, y, family = family, M = 8, cv.flag = cv.flag, verbose = FALSE)
  attr(fit, "long.data") <- list(x = d$features, tm = d$time, id = d$id, y = y)
  fit
}

for (family in c("continuous", "binary", "ordinal", "nominal")) {
  test_that(paste("partial predict matches the full-path computation:", family), {
    fit <- make.family.fit(family, cv.flag = TRUE)
    for (cv in c(FALSE, TRUE)) {
      p <- predict(fit, x = fit$x, tm = fit$time.unique, partial = TRUE, M = 6, use.cv.flag = cv)
      expect_identical(as.q.list(p$mu, family), old.mu(fit, p, p$m.opt))
      if (!cv) {
        expect_identical(
          as.q.list(p$muhat, family),
          old.extrapolate(fit, p$membership, fit$gamma, p$m.opt)
        )
      }
    }
  })
}

test_that("test-response predict keeps the full-path error curve and step choice", {
  for (family in c("ordinal", "nominal")) {
    fit <- make.family.fit(family)
    d <- attr(fit, "long.data")
    p <- predict(fit, x = d$x, tm = d$tm, id = d$id, y = d$y)
    for (m in seq_len(fit$M)) {
      mu.m <- old.mu(fit, p, rep(m, fit$n.q))
      for (q in seq_len(fit$n.q)) {
        expect_identical(
          p$err.rate[[q]][m, ],
          c(
            l1 = boostmtree:::boostmtree.l1.dist(p$y.org[[q]], mu.m[[q]]),
            l2 = boostmtree:::boostmtree.l2.dist(p$y.org[[q]], mu.m[[q]])
          ) / fit$y.sd
        )
      }
    }
    expect_identical(p$mu, old.mu(fit, p, p$m.opt))
  }
})

test_that("every subject on a shared time grid shares one design matrix", {
  fit <- make.family.fit("continuous")
  p <- predict(fit, x = fit$x, tm = fit$time.unique, partial = TRUE)
  expect_identical(p$time.design[[1]], fit$x.tm[seq_len(nrow(fit$x.tm)), , drop = FALSE])
  expect_true(all(vapply(p$time.design, identical, logical(1), p$time.design[[1]])))
})

test_that("partial.plot fills every response component from one prediction per grid value", {
  fit <- make.family.fit("ordinal")
  values <- sort(unique(fit$x$x1))[c(2, 10)]
  for (prob.class in c(FALSE, TRUE)) {
    pp <- partial.plot(
      fit, x.var.names = "x1", x.var.values = list(x1 = values),
      prob.class = prob.class, output = "data"
    )
    time.info <- boostmtree:::boostmtree.effect.time.info(fit)
    for (k in seq_along(values)) {
      new.x <- fit$x
      new.x$x1 <- values[k]
      p <- predict(fit, x = new.x, tm = time.info$observed.time, M = pp$M, partial = TRUE)
      for (q.out in seq_along(pp$response.labels)) {
        response <- if (prob.class) p$prob.class[[q.out]] else p$mu[[q.out]]
        expected <- colMeans(do.call(rbind, lapply(response, function(s) s[time.info$index])))
        expect_identical(unlist(pp$curves[[q.out]]$x1[k, -1], use.names = FALSE), expected)
      }
    }
  }
})
