boostmtree.initialize.lambda <- function(lambda, n.q, univariate, df.time.design, pen.ord, svd.tol) {
  pen.lsq.matrix <- boostmtree.pen.bs.deriv(df.time.design - 1L, pen.ord)
  pen.inv.sqrt.matrix <- NULL
  lambda.est.flag <- FALSE
  if (df.time.design == 1L) {
    lambda <- rep(0, n.q)
  } else if (is.null(lambda)) {
    if (!univariate && df.time.design >= (pen.ord + 2L)) {
      lambda.est.flag <- TRUE
      pen.mix.matrix <- boostmtree.pen.bs(df.time.design - 1L, pen.ord)
      svd.penalty <- svd(pen.mix.matrix)
      d.inverse.sqrt <- 1 / sqrt(svd.penalty$d)
      d.inverse.sqrt[svd.penalty$d < svd.tol] <- 0
      pen.inv.sqrt.matrix <- svd.penalty$v %*% (t(svd.penalty$v) * d.inverse.sqrt)
      lambda <- rep(0, n.q)
    } else {
      if (!univariate) {
        warning(
          "Not enough degrees of freedom to estimate `lambda`; setting it to zero."
        )
      }
      lambda <- rep(0, n.q)
    }
  } else {
    if (length(lambda) == 1L) {
      lambda <- rep(as.numeric(lambda), n.q)
    } else if (length(lambda) != n.q) {
      stop("`lambda` must be NULL, a scalar, or have length equal to `n.q`.")
    } else {
      lambda <- as.numeric(lambda)
    }
    if (any(is.na(lambda)) || any(lambda < 0)) {
      stop("`lambda` must be non-negative when supplied.")
    }
  }
  list(
    lambda = lambda,
    lambda.est.flag = lambda.est.flag,
    pen.lsq.matrix = pen.lsq.matrix,
    pen.inv.sqrt.matrix = pen.inv.sqrt.matrix
  )
}
boostmtree.initialize.rho <- function(rho, n.q, univariate) {
  if (univariate) {
    return(list(
      rho = rep(0, n.q),
      rho.fit.flag = FALSE
    ))
  }
  if (is.null(rho)) {
    return(list(
      rho = rep(0, n.q),
      rho.fit.flag = TRUE
    ))
  }
  if (length(rho) == 1L) {
    rho <- rep(as.numeric(rho), n.q)
  } else if (length(rho) != n.q) {
    stop("`rho` must be NULL, a scalar, or have length equal to `n.q`.")
  } else {
    rho <- as.numeric(rho)
  }
  if (any(is.na(rho)) || any(abs(rho) >= 1)) {
    stop("`rho` must lie strictly inside (-1, 1).")
  }
  list(
    rho = rho,
    rho.fit.flag = FALSE
  )
}
boostmtree.initialize.cv.state <- function(model.info) {
  n.subject <- model.info$n
  n.q <- model.info$n.q
  family <- model.info$family
  mu.cv.list <- lapply(seq_len(n.q), function(q) {
    vector("list", model.info$M)
  })
  if (family == "nominal") {
    l.pred.ref.cv <- lapply(seq_len(n.subject), function(i) {
      rep(log(1 / model.info$q.total), model.info$ni[i])
    })
  } else {
    l.pred.ref.cv <- lapply(seq_len(n.subject), function(i) {
      rep(0, model.info$ni[i])
    })
  }
  l.pred.cv <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n.subject), function(i) {
      l.pred.ref.cv[[i]] + model.info$l.pred.db[[q]][[i]]
    })
  })
  mu.cv <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n.subject), function(i) {
      boostmtree.get.mu(l.pred.cv[[q]][[i]], family = family)
    })
  })
  l.pred.db.i <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n.subject), function(i) {
      lapply(seq_len(n.subject), function(j) {
        rep(0, model.info$ni[j])
      })
    })
  })
  gamma.i.list <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(model.info$M), function(m) {
      vector("list", length = n.subject)
    })
  })
  err.rate <- lapply(seq_len(n.q), function(q) {
    err.rate.matrix <- matrix(NA_real_, model.info$M, 2L)
    colnames(err.rate.matrix) <- c("l1", "l2")
    err.rate.matrix
  })
  if (family == "continuous") {
    y.mean.i <- lapply(seq_len(n.q), function(q) {
      vapply(seq_len(n.subject), function(i) {
        mean(unlist(model.info$y.org[[q]][-i]), na.rm = TRUE)
      }, numeric(1))
    })
    y.sd.i <- lapply(seq_len(n.q), function(q) {
      vapply(seq_len(n.subject), function(i) {
        sd.i <- sd(unlist(model.info$y.org[[q]][-i]), na.rm = TRUE)
        if (sd.i < 1e-6) 1 else sd.i
      }, numeric(1))
    })
  } else {
    y.mean.i <- lapply(seq_len(n.q), function(q) rep(model.info$y.mean, n.subject))
    y.sd.i <- lapply(seq_len(n.q), function(q) rep(model.info$y.sd, n.subject))
  }
  list(
    mu.cv.list = mu.cv.list,
    l.pred.ref.cv = l.pred.ref.cv,
    l.pred.cv = l.pred.cv,
    mu.cv = mu.cv,
    l.pred.db.i = l.pred.db.i,
    gamma.i.list = gamma.i.list,
    err.rate = err.rate,
    m.opt = rep(NA_integer_, n.q),
    rmse = rep(NA_real_, n.q),
    y.mean.i = y.mean.i,
    y.sd.i = y.sd.i
  )
}
boostmtree.make.variance.matrices <- function(ni, rho, phi) {
  v.matrix <- lapply(seq_along(ni), function(i) {
    variance.matrix <- matrix(rho * phi, ni[i], ni[i])
    diag(variance.matrix) <- phi
    variance.matrix
  })
  inv.v.matrix <- lapply(seq_along(v.matrix), function(i) {
    out <- tryCatch(qr.solve(v.matrix[[i]]), error = function(ex) NULL)
    if (is.null(out)) {
      diag(1 / max(phi, .Machine$double.eps), nrow(v.matrix[[i]]))
    } else {
      out
    }
  })
  list(v.matrix = v.matrix, inv.v.matrix = inv.v.matrix)
}
boostmtree.compute.working.gradient <- function(y.current, mu.current, time.design, inv.v.matrix, family, mod.grad) {
  df.time.design <- ncol(time.design[[1]])
  n.subject <- length(y.current)
  h.mu <- lapply(seq_len(n.subject), function(i) {
    boostmtree.transform.h(mu.current[[i]], family = family)
  })
  if (!mod.grad) {
    gradient.list <- lapply(seq_len(n.subject), function(i) {
      t(time.design[[i]]) %*% h.mu[[i]] %*% inv.v.matrix[[i]] %*% (y.current[[i]] - mu.current[[i]])
    })
  } else {
    gradient.list <- lapply(seq_len(n.subject), function(i) {
      t(time.design[[i]]) %*% h.mu[[i]] %*% (y.current[[i]] - mu.current[[i]])
    })
  }
  t(matrix(unlist(gradient.list), nrow = df.time.design))
}
boostmtree.fit.gamma.cluster <- function(subject.index, l.pred.current, y.current, time.design, inv.v.matrix, lambda, pen.lsq.matrix, family, nr.iter) {
  df.time.design <- ncol(time.design[[1]])
  if (length(subject.index) == 0L) {
    return(rep(0, df.time.design))
  }
  gamma.update <- rep(0, df.time.design)
  for (iter in seq_len(nr.iter)) {
    mu.update <- lapply(subject.index, function(i) {
      boostmtree.get.mu(
        linear.predictor = l.pred.current[[i]] + c(time.design[[i]] %*% gamma.update),
        family = family
      )
    })
    cal.d <- lapply(seq_along(subject.index), function(j) {
      boostmtree.transform.h(mu.update[[j]], family = family) %*% time.design[[subject.index[j]]]
    })
    hessian.temp <- Reduce("+", lapply(seq_along(subject.index), function(j) {
      t(cal.d[[j]]) %*% inv.v.matrix[[subject.index[j]]] %*% cal.d[[j]]
    }))
    score.temp <- Reduce("+", lapply(seq_along(subject.index), function(j) {
      t(cal.d[[j]]) %*% inv.v.matrix[[subject.index[j]]] %*% (
        y.current[[subject.index[j]]] - mu.update[[j]]
      )
    }))
    hessian <- hessian.temp + lambda * pen.lsq.matrix
    score <- score.temp - lambda * (pen.lsq.matrix %*% gamma.update)
    qr.object <- tryCatch(qr.solve(hessian, score), error = function(ex) NULL)
    if (is.null(qr.object)) {
      qr.object <- rep(0, df.time.design)
    }
    gamma.update <- gamma.update + qr.object
  }
  gamma.update
}
boostmtree.make.gamma.step <- function(gamma.coefficients, membership.original) {
  list(
    node.label = sort(unique(membership.original)),
    coefficients = matrix(
      unlist(gamma.coefficients),
      ncol = length(gamma.coefficients[[1]]),
      byrow = TRUE
    )
  )
}
boostmtree.expand.gamma.by.subject <- function(gamma.step, membership.original) {
  coefficient.count <- ncol(gamma.step$coefficients)
  matrix(
    unlist(lapply(seq_along(membership.original), function(i) {
      row.index <- match(
        as.character(membership.original[i]),
        as.character(gamma.step$node.label)
      )
      gamma.step$coefficients[row.index, ]
    })),
    ncol = coefficient.count,
    byrow = TRUE
  )
}
boostmtree.estimate.lambda <- function(membership, node.count, y.current, mu.current, l.pred.current, time.design, rho, phi, family, pen.inv.sqrt.matrix, lambda.initial, lambda.iter, lambda.max) {
  transformed.data <- lapply(seq_along(y.current), function(i) {
    n.i <- length(y.current[[i]])
    if (n.i > 1L) {
      c.i <- boostmtree.rho.inv.sqrt(n.i, rho)
      r.inv.sqrt <- (diag(1, n.i) - matrix(c.i, n.i, n.i)) / sqrt(1 - rho)
      v.inv.sqrt <- phi^(-1 / 2) * r.inv.sqrt
    } else {
      r.inv.sqrt <- cbind(1)
      v.inv.sqrt <- phi^(-1 / 2) * r.inv.sqrt
    }
    y.new <- v.inv.sqrt %*% (y.current[[i]] - mu.current[[i]])
    mu.lambda <- boostmtree.get.mu.lambda(l.pred.current[[i]], family = family)
    lambda.design <- boostmtree.transform.h(mu.lambda, family = family) %*% time.design[[i]]
    x.new <- v.inv.sqrt %*% lambda.design[, 1, drop = FALSE]
    z.new <- v.inv.sqrt %*% lambda.design[, -1, drop = FALSE] %*% pen.inv.sqrt.matrix
    list(y.new = y.new, x.new = x.new, z.new = z.new)
  })
  lambda.hat <- lambda.initial
  for (iteration in seq_len(lambda.iter)) {
    blup.object <- boostmtree.blup.solve(
      transformed.data = transformed.data,
      membership = membership,
      sigma = lambda.hat,
      node.count = node.count
    )
    lambda.object <- lapply(seq_len(node.count), function(node.id) {
      node.subject <- membership == node.id
      z.matrix <- do.call(rbind, lapply(which(node.subject), function(j) transformed.data[[j]]$z.new))
      x.matrix <- do.call(rbind, lapply(which(node.subject), function(j) transformed.data[[j]]$x.new))
      y.vector <- unlist(lapply(which(node.subject), function(j) transformed.data[[j]]$y.new))
      zz <- t(z.matrix) %*% z.matrix
      rss <- (y.vector - x.matrix %*% c(blup.object[[node.id]]$fixed.effect))^2
      robust.point <- rss <= quantile(rss, 0.99, na.rm = TRUE)
      rss <- sum(rss[robust.point], na.rm = TRUE)
      residual <- (
        y.vector -
          x.matrix %*% c(blup.object[[node.id]]$fixed.effect) -
          z.matrix %*% c(blup.object[[node.id]]$random.effect)
      )^2
      residual <- residual[robust.point]
      list(
        trace.z = sum(diag(zz)),
        rss = rss,
        residual = residual
      )
    })
    numerator <- sum(vapply(lambda.object, function(obj) obj$trace.z, numeric(1)), na.rm = TRUE)
    denominator <- sum(vapply(lambda.object, function(obj) obj$rss, numeric(1)), na.rm = TRUE)
    n.residual <- sum(unlist(lapply(lambda.object, function(obj) obj$residual)), na.rm = TRUE)
    if (!is.na(denominator) && denominator > (0.99 * n.residual)) {
      lambda.hat <- numerator / (denominator - 0.99 * n.residual)
    } else {
      lambda.hat <- min(lambda.hat, lambda.max)
    }
    lambda.hat <- min(lambda.hat, lambda.max)
  }
  lambda.hat
}
boostmtree.update.cv.step <- function(cv.state.q, rfsrc.obj, membership, membership.original, time.design, nu.vec, lambda, y.current, inv.v.matrix, family, nr.iter, pen.lsq.matrix, previous.ordinal = NULL) {
  n.subject <- length(cv.state.q)
  node.count <- length(unique(membership))
  df.time.design <- ncol(time.design[[1]])
  gamma.at.m <- vector("list", n.subject)
  updated <- cv.state.q
  oob <- which(rfsrc.obj$inbag == 0)
  if (length(oob) > 0L) {
    for (i in oob) {
      l.pred.ij <- cv.state.q[[i]]
      gamma.coefficients <- lapply(seq_len(node.count), function(node.id) {
        subject.index <- setdiff(which(membership == node.id), i)
        if (length(subject.index) == 0L) {
          rep(0, df.time.design)
        } else {
          boostmtree.fit.gamma.cluster(
            subject.index = subject.index,
            l.pred.current = l.pred.ij,
            y.current = y.current,
            time.design = time.design,
            inv.v.matrix = inv.v.matrix,
            lambda = lambda,
            pen.lsq.matrix = pen.lsq.matrix,
            family = family,
            nr.iter = nr.iter
          )
        }
      })
      gamma.step.i <- boostmtree.make.gamma.step(gamma.coefficients, membership.original)
      gamma.at.m[[i]] <- gamma.step.i
      updated[[i]] <- lapply(seq_len(n.subject), function(j) {
        row.index <- match(
          as.character(membership.original[j]),
          as.character(gamma.step.i$node.label)
        )
        l.pred.ij[[j]] + c(time.design[[j]] %*% (gamma.step.i$coefficients[row.index, ] * nu.vec))
      })
    }
  }
  if (family == "ordinal" && !is.null(previous.ordinal)) {
    updated <- lapply(seq_len(n.subject), function(i) {
      lapply(seq_len(n.subject), function(j) {
        pmax(updated[[i]][[j]], previous.ordinal[[i]][[j]])
      })
    })
  }
  list(
    l.pred.db.i.q = updated,
    gamma.i.at.m = gamma.at.m
  )
}
boostmtree.update.gls.parameters <- function(x.long, time, id, y.current, mu.current, rho, phi) {
  residual.data <- data.frame(
    y = unlist(lapply(seq_along(y.current), function(i) y.current[[i]] - mu.current[[i]])),
    x.long,
    tm = time,
    id = id
  )
  gls.object <- tryCatch(
    gls(y ~ ., data = residual.data, correlation = corCompSymm(form = ~ 1 | id)),
    error = function(ex) NULL
  )
  if (is.null(gls.object)) {
    gls.object <- tryCatch(
      gls(y ~ 1, data = residual.data, correlation = corCompSymm(form = ~ 1 | id)),
      error = function(ex) NULL
    )
  }
  if (is.null(gls.object)) {
    return(list(
      rho = rho,
      phi = phi
    ))
  }
  cor.object <- gls.object$modelStruct$corStruct
  if (is.null(cor.object)) {
    cor.object <- gls.object$modelStruct$corStruc
  }
  rho.temp <- as.numeric(coef(cor.object, unconstrained = FALSE))
  list(
    rho = max(min(0.999, rho.temp, na.rm = TRUE), -0.999),
    phi = gls.object$sigma^2
  )
}
boostmtree.finalize.fit <- function(model.info, fit.info, cv.state = NULL) {
  if (!isTRUE(model.info$cv.flag)) {
    mu.final <- lapply(seq_len(model.info$n.q), function(q) {
      lapply(seq_len(model.info$n), function(i) {
        boostmtree.get.mu(
          linear.predictor = fit.info$l.pred[[q]][[i]] * model.info$y.sd + model.info$y.mean,
          family = model.info$family
        )
      })
    })
    return(list(
      mu = mu.final,
      err.rate = NULL,
      rmse = NULL,
      m.opt = NULL,
      gamma.i.list = NULL
    ))
  }
  for (q in seq_len(model.info$n.q)) {
    diff.err <- abs(
      cv.state$err.rate[[q]][, "l2"] -
        min(cv.state$err.rate[[q]][, "l2"], na.rm = TRUE)
    )
    diff.err[is.na(diff.err)] <- Inf
    threshold <- model.info$y.sd * model.info$eps
    acceptable <- which(diff.err < threshold)
    if (length(acceptable) > 0L) {
      cv.state$m.opt[q] <- min(acceptable)
    } else {
      cv.state$m.opt[q] <- model.info$M
    }
    cv.state$rmse[q] <- cv.state$err.rate[[q]][cv.state$m.opt[q], "l2"]
    fit.info$mu[[q]] <- cv.state$mu.cv.list[[q]][[cv.state$m.opt[q]]]
  }
  list(
    mu = fit.info$mu,
    err.rate = cv.state$err.rate,
    rmse = cv.state$rmse,
    m.opt = cv.state$m.opt,
    gamma.i.list = cv.state$gamma.i.list
  )
}
boostmtree.fit.tree <- function(model.info) {
  if (model.info$family == "continuous") {
    model.info$nr.iter <- 1L
  }
  if (length(model.info$nu) == 1L) {
    model.info$nu <- rep(model.info$nu, 2L)
  }
  if (length(model.info$nu) != 2L || any(!(0 < model.info$nu & model.info$nu <= 1))) {
    stop("`nu` must have length 1 or 2, with all values in (0, 1].")
  }
  model.info$mtry <- if (is.null(model.info$mtry)) {
    ncol(model.info$x.subject)
  } else {
    as.integer(model.info$mtry)
  }
  if (length(model.info$mtry) != 1L || is.na(model.info$mtry) || model.info$mtry < 1L) {
    stop("`mtry` must be NULL or a single positive integer.")
  }
  model.info$k <- as.integer(model.info$k)
  if (length(model.info$k) != 1L || is.na(model.info$k) || model.info$k < 1L) {
    stop("`k` must be a single positive integer.")
  }
  model.info$lambda.iter <- as.integer(model.info$lambda.iter)
  model.info$nr.iter <- as.integer(model.info$nr.iter)
  if (model.info$lambda.iter < 1L || model.info$nr.iter < 1L) {
    stop("`lambda.iter` and `nr.iter` must be positive integers.")
  }
  lambda.info <- boostmtree.initialize.lambda(
    lambda = model.info$lambda,
    n.q = model.info$n.q,
    univariate = model.info$univariate,
    df.time.design = model.info$df.time.design,
    pen.ord = model.info$pen.ord,
    svd.tol = model.info$svd.tol
  )
  rho.info <- boostmtree.initialize.rho(
    rho = model.info$rho,
    n.q = model.info$n.q,
    univariate = model.info$univariate
  )
  lambda <- lambda.info$lambda
  rho <- rho.info$rho
  rho.fit.flag <- rho.info$rho.fit.flag
  lambda.est.flag <- lambda.info$lambda.est.flag
  pen.lsq.matrix <- lambda.info$pen.lsq.matrix
  pen.inv.sqrt.matrix <- lambda.info$pen.inv.sqrt.matrix
  model.info$nu.vec <- c(
    model.info$nu[1],
    rep(model.info$nu[2], model.info$df.time.design - 1L)
  )
  sigma <- phi <- rep(1, model.info$n.q)
  if (!lambda.est.flag) {
    sigma <- vapply(seq_len(model.info$n.q), function(q) {
      boostmtree.sigma.robust(lambda[q], rho[q])
    }, numeric(1))
  }
  if (model.info$family == "nominal") {
    l.pred.ref <- lapply(seq_len(model.info$n), function(i) {
      rep(log(1 / model.info$q.total), model.info$ni[i])
    })
  } else {
    l.pred.ref <- lapply(seq_len(model.info$n), function(i) {
      rep(0, model.info$ni[i])
    })
  }
  l.pred.db <- lapply(seq_len(model.info$n.q), function(q) {
    lapply(seq_len(model.info$n), function(i) rep(0, model.info$ni[i]))
  })
  model.info$l.pred.db <- l.pred.db
  y.by.q <- lapply(seq_len(model.info$n.q), function(q) {
    boostmtree.split.by.subject(
      model.info$y.by.q.vector[[q]],
      id = model.info$id,
      id.unique = model.info$id.unique
    )
  })
  y.standardized <- lapply(seq_len(model.info$n.q), function(q) {
    lapply(seq_len(model.info$n), function(i) {
      (y.by.q[[q]][[i]] - model.info$y.mean) / model.info$y.sd
    })
  })
  base.learner <- lapply(seq_len(model.info$n.q), function(q) vector("list", model.info$M))
  membership.list <- lapply(seq_len(model.info$n.q), function(q) vector("list", model.info$M))
  gamma.list <- lapply(seq_len(model.info$n.q), function(q) vector("list", model.info$M))
  oob.subject.count <- matrix(0L, nrow = model.info$M, ncol = model.info$n.q)
  if (!model.info$univariate) {
    lambda.matrix <- matrix(NA_real_, nrow = model.info$M, ncol = model.info$n.q)
    phi.matrix <- matrix(NA_real_, nrow = model.info$M, ncol = model.info$n.q)
    rho.matrix <- matrix(NA_real_, nrow = model.info$M, ncol = model.info$n.q)
  } else {
    lambda.matrix <- NULL
    phi.matrix <- NULL
    rho.matrix <- NULL
  }
  cv.state <- NULL
  model.info$cv.flag <- isTRUE(model.info$cv.flag)
  cv.lambda.flag <- model.info$cv.flag && isTRUE(model.info$control$cv.lambda) && lambda.est.flag
  cv.rho.flag <- model.info$cv.flag && isTRUE(model.info$control$cv.rho) && rho.fit.flag
  if (model.info$cv.flag) {
    cv.state <- boostmtree.initialize.cv.state(model.info = model.info)
  }
  y.names <- paste0("y", seq_len(model.info$df.time.design))
  rfsrc.formula <- as.formula(
    paste("Multivar(", paste(y.names, collapse = ","), ") ~ .", sep = "")
  )
  sample.matrix <- boostmtree.resolve.sample.matrix(
    control = model.info$control,
    n = model.info$n,
    M = model.info$M
  )
  l.pred <- lapply(seq_len(model.info$n.q), function(q) {
    lapply(seq_len(model.info$n), function(i) l.pred.ref[[i]] + l.pred.db[[q]][[i]])
  })
  mu <- lapply(seq_len(model.info$n.q), function(q) {
    lapply(seq_len(model.info$n), function(i) {
      boostmtree.get.mu(l.pred[[q]][[i]], family = model.info$family)
    })
  })
  lambda.initial <- model.info$y.sd^2
  nodesize <- max(1L, round(model.info$n / (2 * model.info$k)))
  progress.bar <- NULL
  if (isTRUE(model.info$verbose)) {
    progress.bar <- txtProgressBar(min = 0, max = model.info$M, style = 3)
    on.exit(close(progress.bar), add = TRUE)
  }
  for (m in seq_len(model.info$M)) {
    if (isTRUE(model.info$verbose)) {
      setTxtProgressBar(progress.bar, m)
      if (m == model.info$M) {
        cat("\n")
      }
    }
    for (q in seq_len(model.info$n.q)) {
      v.matrix <- lapply(seq_len(model.info$n), function(i) {
        variance.matrix <- matrix(rho[q] * phi[q], model.info$ni[i], model.info$ni[i])
        diag(variance.matrix) <- phi[q]
        variance.matrix
      })
      inv.v.matrix <- lapply(seq_len(model.info$n), function(i) {
        out <- tryCatch(qr.solve(v.matrix[[i]]), error = function(ex) NULL)
        if (is.null(out)) {
          diag(1 / max(phi[q], .Machine$double.eps), nrow(v.matrix[[i]]))
        } else {
          out
        }
      })
      h.mu <- lapply(seq_len(model.info$n), function(i) {
        boostmtree.transform.h(mu[[q]][[i]], family = model.info$family)
      })
      if (!model.info$mod.grad) {
        gm.mod <- t(matrix(unlist(lapply(seq_len(model.info$n), function(i) {
          t(model.info$time.design[[i]]) %*% h.mu[[i]] %*% inv.v.matrix[[i]] %*% (
            y.standardized[[q]][[i]] - mu[[q]][[i]]
          )
        })), nrow = model.info$df.time.design))
      } else {
        gm.mod <- t(matrix(unlist(lapply(seq_len(model.info$n), function(i) {
          t(model.info$time.design[[i]]) %*% h.mu[[i]] %*% (
            y.standardized[[q]][[i]] - mu[[q]][[i]]
          )
        })), nrow = model.info$df.time.design))
      }
      incoming.data <- cbind(gm.mod, model.info$x.subject)
      incoming.data <- as.data.frame(incoming.data, check.names = FALSE)
      names(incoming.data) <- c(y.names, names(model.info$x.subject))
      rfsrc.object <- rfsrc(
        rfsrc.formula,
        data = incoming.data,
        ntree = 1,
        mtry = model.info$mtry,
        nodesize = nodesize,
        nsplit = model.info$control$nsplit,
        importance = "none",
        bootstrap = model.info$control$bootstrap,
        samptype = model.info$control$samptype,
        samp = if (identical(model.info$control$bootstrap, "by.user")) sample.matrix[, m, drop = FALSE] else NULL,
        xvar.wt = model.info$control$xvar.wt,
        case.wt = model.info$control$case.wt,
        membership = TRUE,
        na.action = model.info$na.action,
        nimpute = 1,
        seed = model.info$control$seed
      )
      base.learner[[q]][[m]] <- rfsrc.object
      oob.subject.count[m, q] <- sum(rfsrc.object$inbag == 0L)
      prediction.object <- predict.rfsrc(
        rfsrc.object,
        membership = TRUE,
        ptn.count = model.info$k,
        importance = "none"
      )
      membership.original <- c(prediction.object$ptn.membership)
      membership.list[[q]][[m]] <- membership.original
      membership <- as.numeric(factor(membership.original))
      node.count <- length(unique(membership))
      if (lambda.est.flag) {
        mu.for.lambda <- if (cv.lambda.flag) cv.state$mu.cv[[q]] else mu[[q]]
        lambda[q] <- boostmtree.estimate.lambda(
          membership = membership,
          node.count = node.count,
          y.current = y.standardized[[q]],
          mu.current = mu.for.lambda,
          l.pred.current = l.pred[[q]],
          time.design = model.info$time.design,
          rho = rho[q],
          phi = phi[q],
          family = model.info$family,
          pen.inv.sqrt.matrix = pen.inv.sqrt.matrix,
          lambda.initial = lambda.initial,
          lambda.iter = model.info$lambda.iter,
          lambda.max = model.info$lambda.max
        )
        sigma[q] <- boostmtree.sigma.robust(lambda[q], rho[q])
      }
      gamma.coefficients <- lapply(seq_len(node.count), function(node.id) {
        node.subject <- membership == node.id
        if (sum(node.subject) > 0) {
          subject.index <- which(node.subject)
          subject.sequence <- seq_along(subject.index)
          gamma.update <- rep(0, model.info$df.time.design)
          for (iter in seq_len(model.info$nr.iter)) {
            mu.update <- lapply(subject.index, function(i.subject) {
              l.pred.gamma <- l.pred[[q]][[i.subject]] + c(model.info$time.design[[i.subject]] %*% gamma.update)
              boostmtree.get.mu(l.pred.gamma, family = model.info$family)
            })
            cal.d <- lapply(subject.sequence, function(j) {
              h.mu.update <- boostmtree.transform.h(mu.update[[j]], family = model.info$family)
              h.mu.update %*% model.info$time.design[[subject.index[j]]]
            })
            hessian.temp <- Reduce("+", lapply(subject.sequence, function(j) {
              t(cal.d[[j]]) %*% inv.v.matrix[[subject.index[j]]] %*% cal.d[[j]]
            }))
            score.temp <- Reduce("+", lapply(subject.sequence, function(j) {
              t(cal.d[[j]]) %*% inv.v.matrix[[subject.index[j]]] %*% (
                y.standardized[[q]][[subject.index[j]]] - mu.update[[j]]
              )
            }))
            hessian <- hessian.temp + lambda[q] * pen.lsq.matrix
            score <- score.temp - lambda[q] * (pen.lsq.matrix %*% gamma.update)
            qr.object <- tryCatch(qr.solve(hessian, score), error = function(ex) NULL)
            if (is.null(qr.object)) {
              qr.object <- rep(0, model.info$df.time.design)
            }
            gamma.update <- gamma.update + qr.object
          }
          gamma.update
        } else {
          rep(0, model.info$df.time.design)
        }
      })
      gamma.matrix <- matrix(0, nrow = node.count, ncol = model.info$df.time.design + 1L)
      gamma.matrix[, 1L] <- sort(unique(membership.original))
      gamma.matrix[, 2:(model.info$df.time.design + 1L)] <- matrix(
        unlist(gamma.coefficients),
        ncol = model.info$df.time.design,
        byrow = TRUE
      )
      gamma.list[[q]][[m]] <- gamma.matrix
      gamma.by.subject <- t(
        matrix(
          unlist(lapply(seq_len(model.info$n), function(i) {
            gamma.coefficients[[membership[i]]]
          })),
          nrow = model.info$df.time.design
        ) * model.info$nu.vec
      )
      l.pred.db[[q]] <- lapply(seq_len(model.info$n), function(i) {
        updated <- l.pred.db[[q]][[i]] + c(model.info$time.design[[i]] %*% gamma.by.subject[i, ])
        if (model.info$family == "ordinal" && q > 1L) {
          pmax(updated, l.pred.db[[q - 1L]][[i]])
        } else {
          updated
        }
      })
      if (model.info$cv.flag) {
        cv.update <- boostmtree.update.cv.step(
          cv.state.q = cv.state$l.pred.db.i[[q]],
          rfsrc.obj = rfsrc.object,
          membership = membership,
          membership.original = membership.original,
          time.design = model.info$time.design,
          nu.vec = model.info$nu.vec,
          lambda = lambda[q],
          y.current = y.standardized[[q]],
          inv.v.matrix = inv.v.matrix,
          family = model.info$family,
          nr.iter = model.info$nr.iter,
          pen.lsq.matrix = pen.lsq.matrix,
          previous.ordinal = if (model.info$family == "ordinal" && q > 1L) cv.state$l.pred.db.i[[q - 1L]] else NULL
        )
        cv.state$l.pred.db.i[[q]] <- cv.update$l.pred.db.i.q
        cv.state$gamma.i.list[[q]][[m]] <- cv.update$gamma.i.at.m
      }
    }
    if (model.info$cv.flag) {
      if (model.info$family == "nominal") {
        cv.state$l.pred.ref.cv <- lapply(seq_len(model.info$n), function(i) {
          log((1 + Reduce("+", lapply(seq_len(model.info$n.q), function(q) {
            exp(cv.state$l.pred.db.i[[q]][[i]][[i]])
          })))^(-1))
        })
      }
      for (q in seq_len(model.info$n.q)) {
        cv.state$l.pred.cv[[q]] <- lapply(seq_len(model.info$n), function(i) {
          cv.state$l.pred.ref.cv[[i]] + cv.state$l.pred.db.i[[q]][[i]][[i]]
        })
        cv.state$mu.cv[[q]] <- lapply(seq_len(model.info$n), function(i) {
          boostmtree.get.mu(cv.state$l.pred.cv[[q]][[i]], family = model.info$family)
        })
        l.pred.cv.original.scale <- lapply(seq_len(model.info$n), function(i) {
          cv.state$l.pred.cv[[q]][[i]] * cv.state$y.sd.i[[q]][i] + cv.state$y.mean.i[[q]][i]
        })
        mu.cv.original.scale <- lapply(seq_len(model.info$n), function(i) {
          boostmtree.get.mu(l.pred.cv.original.scale[[i]], family = model.info$family)
        })
        cv.state$mu.cv.list[[q]][[m]] <- mu.cv.original.scale
        cv.state$err.rate[[q]][m, ] <- c(
          boostmtree.l1.dist(model.info$y.org[[q]], mu.cv.original.scale),
          boostmtree.l2.dist(model.info$y.org[[q]], mu.cv.original.scale)
        )
      }
    } else {
      if (model.info$family == "nominal") {
        l.pred.ref <- lapply(seq_len(model.info$n), function(i) {
          log((1 + Reduce("+", lapply(seq_len(model.info$n.q), function(q) {
            exp(l.pred.db[[q]][[i]])
          })))^(-1))
        })
      }
      for (q in seq_len(model.info$n.q)) {
        l.pred[[q]] <- lapply(seq_len(model.info$n), function(i) {
          l.pred.ref[[i]] + l.pred.db[[q]][[i]]
        })
        mu[[q]] <- lapply(seq_len(model.info$n), function(i) {
          boostmtree.get.mu(l.pred[[q]][[i]], family = model.info$family)
        })
      }
    }
    for (q in seq_len(model.info$n.q)) {
      if (!model.info$univariate && rho.fit.flag) {
        gls.current <- if (cv.rho.flag) cv.state$mu.cv[[q]] else mu[[q]]
        gls.update <- boostmtree.update.gls.parameters(
          x.long = model.info$x.long,
          time = model.info$time,
          id = model.info$id,
          y.current = y.standardized[[q]],
          mu.current = gls.current,
          rho = rho[q],
          phi = phi[q]
        )
        rho[q] <- gls.update$rho
        phi[q] <- gls.update$phi
      }
      if (!model.info$univariate) {
        phi.matrix[m, q] <- phi[q] * model.info$y.sd^2
        rho.matrix[m, q] <- rho[q]
        sigma[q] <- boostmtree.sigma.robust(lambda[q], rho[q])
        lambda.matrix[m, q] <- lambda[q]
      }
    }
  }
  if (isTRUE(model.info$cv.flag) && any(oob.subject.count == 0L)) {
    warning(
      "One or more boosting iterations had no out-of-bag subjects. ",
      "Cross-validation and grow-object variable importance may be unstable for those iterations."
    )
  }
  final.fit <- boostmtree.finalize.fit(
    model.info = model.info,
    fit.info = list(
      l.pred = l.pred,
      mu = mu
    ),
    cv.state = cv.state
  )
  list(
    y.org = model.info$y.org,
    mu = final.fit$mu,
    lambda = lambda.matrix,
    phi = phi.matrix,
    rho = rho.matrix,
    base.learner = base.learner,
    membership = membership.list,
    gamma = gamma.list,
    err.rate = final.fit$err.rate,
    rmse = final.fit$rmse,
    m.opt = final.fit$m.opt,
    gamma.i.list = final.fit$gamma.i.list,
    oob.subject.count = oob.subject.count,
    oob.available = all(oob.subject.count > 0L)
  )
}
