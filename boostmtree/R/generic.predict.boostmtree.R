boostmtree.predict.approx.match <- function(x, y) {
  vapply(seq_along(x), function(i) {
    which.min(abs(x[i] - y))
  }, integer(1))
}
boostmtree.predict.q.set.index <- function(object) {
  if (identical(object$family, "continuous")) {
    return(NA_integer_)
  }
  match(object$q.set, object$y.levels)
}
boostmtree.predict.validate.x <- function(x, training.names) {
  x <- as.data.frame(x)
  if (is.null(colnames(x))) {
    if (ncol(x) != length(training.names)) {
      stop(
        "`x` must contain ", length(training.names),
        " columns to match the fitted model."
      )
    }
    colnames(x) <- training.names
  } else {
    missing.columns <- setdiff(training.names, colnames(x))
    if (length(missing.columns) > 0L) {
      stop(
        "`x` is missing required covariates: ",
        paste(missing.columns, collapse = ", ")
      )
    }
    x <- x[, training.names, drop = FALSE]
  }
  x
}
boostmtree.predict.prepare.long.data <- function(x, tm, id, y = NULL, na.action, training.names) {
  x <- boostmtree.predict.validate.x(x, training.names = training.names)
  if (length(tm) != nrow(x)) {
    stop("`tm` must have length equal to `nrow(x)`.")
  }
  if (length(id) != nrow(x)) {
    stop("`id` must have length equal to `nrow(x)`.")
  }
  if (!is.null(y) && length(y) != nrow(x)) {
    stop("`y` must have length equal to `nrow(x)` when supplied.")
  }
  id.unique <- sort(unique(id))
  row.order <- order(match(id, id.unique), seq_along(id))
  x.long <- x[row.order, , drop = FALSE]
  tm <- tm[row.order]
  id <- id[row.order]
  if (!is.null(y)) {
    y <- y[row.order]
  }
  subject.start <- match(id.unique, id)
  x.subject.raw <- x.long[subject.start, , drop = FALSE]
  missing.info <- boostmtree.remove.missing.subjects(
    x.subject.raw,
    na.action = na.action
  )
  kept.ids <- id.unique[missing.info$keep.subject]
  keep.long <- id %in% kept.ids
  x.long <- x.long[keep.long, missing.info$keep.columns, drop = FALSE]
  tm <- tm[keep.long]
  id <- id[keep.long]
  if (!is.null(y)) {
    y <- y[keep.long]
  }
  id.unique <- kept.ids
  subject.start <- match(id.unique, id)
  x.subject <- x.long[subject.start, , drop = FALSE]
  list(
    x = x.subject,
    time = boostmtree.split.by.subject(tm, id, id.unique),
    id = id,
    id.unique = id.unique,
    y = y
  )
}
boostmtree.predict.encode.response <- function(y, object, id, id.unique) {
  family <- object$family
  if (identical(family, "continuous")) {
    if (!is.numeric(y)) {
      stop("`y` must be numeric for the continuous family.")
    }
    y.raw <- boostmtree.split.by.subject(y, id, id.unique)
    return(list(
      y = y.raw,
      y.org = list(y.raw)
    ))
  }
  y.levels <- object$y.levels
  q.set.index <- boostmtree.predict.q.set.index(object)
  y.value <- if (is.factor(y)) {
    as.character(droplevels(y))
  } else if (is.logical(y)) {
    as.character(y)
  } else {
    as.character(y)
  }
  y.code <- match(y.value, as.character(y.levels))
  if (anyNA(y.code)) {
    stop("Observed response values in `y` do not match the response levels used to fit `object`.")
  }
  y.raw <- boostmtree.split.by.subject(y, id, id.unique)
  y.org <- lapply(seq_len(object$n.q), function(q) {
    y.q <- if (family %in% c("binary", "nominal")) {
      as.integer(y.code == q.set.index[q])
    } else {
      as.integer(y.code <= q.set.index[q])
    }
    boostmtree.split.by.subject(y.q, id, id.unique)
  })
  list(
    y = y.raw,
    y.org = y.org
  )
}
boostmtree.predict.normalize.gamma.step <- function(gamma.step) {
  if (is.list(gamma.step) && !is.null(gamma.step$node.label) && !is.null(gamma.step$coefficients)) {
    coefficient.matrix <- gamma.step$coefficients
    if (is.null(dim(coefficient.matrix))) {
      coefficient.matrix <- matrix(coefficient.matrix, nrow = 1L)
    }
    return(list(
      node.label = gamma.step$node.label,
      coefficients = coefficient.matrix
    ))
  }
  if (is.matrix(gamma.step)) {
    if (ncol(gamma.step) < 2L) {
      stop("Stored terminal-node coefficients must include a node label column and at least one coefficient column.")
    }
    return(list(
      node.label = gamma.step[, 1, drop = TRUE],
      coefficients = gamma.step[, -1, drop = FALSE]
    ))
  }
  if (is.atomic(gamma.step) && !is.null(gamma.step)) {
    if (length(gamma.step) < 2L) {
      stop("Stored terminal-node coefficients must include a node label and at least one coefficient.")
    }
    return(list(
      node.label = gamma.step[1],
      coefficients = matrix(gamma.step[-1], nrow = 1L)
    ))
  }
  stop("Unsupported terminal-node coefficient format in `gamma`.")
}
boostmtree.predict.match.gamma <- function(gamma.step, membership.label) {
  gamma.step <- boostmtree.predict.normalize.gamma.step(gamma.step)
  row.index <- match(
    as.character(membership.label),
    as.character(gamma.step$node.label)
  )
  if (anyNA(row.index)) {
    stop("Unable to match predicted node memberships to stored terminal-node coefficients.")
  }
  gamma.step$coefficients[row.index, , drop = FALSE]
}
boostmtree.predict.time.design <- function(time.by.subject, object) {
  time.index <- lapply(time.by.subject, function(time.i) {
    boostmtree.predict.approx.match(time.i, object$time.unique)
  })
  lapply(time.index, function(index.i) {
    object$x.tm[index.i, , drop = FALSE]
  })
}
boostmtree.predict.reference.path <- function(l.pred.db.list, family, q.total, n, ni, m.max) {
  if (identical(family, "nominal")) {
    lapply(seq_len(m.max), function(m) {
      lapply(seq_len(n), function(i) {
        log((1 + Reduce("+", lapply(seq_along(l.pred.db.list), function(q) {
          exp(l.pred.db.list[[q]][[m]][[i]])
        })))^(-1))
      })
    })
  } else {
    lapply(seq_len(m.max), function(m) {
      lapply(seq_len(n), function(i) rep(0, ni[i]))
    })
  }
}
boostmtree.predict.select.m <- function(object, err.rate, M, eps) {
  n.q <- object$n.q
  if (!is.null(M)) {
    return(rep(as.integer(M), n.q))
  }
  if (!is.null(err.rate)) {
    m.opt <- rep(NA_integer_, n.q)
    for (q in seq_len(n.q)) {
      diff.err <- abs(err.rate[[q]][, "l2"] - min(err.rate[[q]][, "l2"], na.rm = TRUE))
      diff.err[is.na(diff.err)] <- Inf
      acceptable <- which(diff.err < object$y.sd * eps)
      if (length(acceptable) > 0L) {
        m.opt[q] <- min(acceptable)
      } else {
        m.opt[q] <- nrow(err.rate[[q]])
      }
    }
  } else if (!is.null(object$m.opt)) {
    if (length(object$m.opt) == 1L) {
      m.opt <- rep(as.integer(object$m.opt), n.q)
    } else {
      m.opt <- as.integer(object$m.opt)
    }
  } else {
    m.opt <- rep(object$M, n.q)
  }
  if (identical(object$family, "ordinal")) {
    m.opt <- rep(max(m.opt, na.rm = TRUE), n.q)
  }
  m.opt
}
boostmtree.predict.flatten.single.response <- function(x, family) {
  if (is.null(x)) {
    return(NULL)
  }
  if (family %in% c("nominal", "ordinal")) {
    return(x)
  }
  x[[1L]]
}
boostmtree.predict.flatten.error <- function(x, family) {
  if (is.null(x)) {
    return(NULL)
  }
  if (family %in% c("nominal", "ordinal")) {
    return(x)
  }
  x[[1L]]
}
boostmtree.predict.extrapolate <- function(object, membership, gamma, m.opt) {
  n <- if (!is.null(membership[[1]][[1]])) length(membership[[1]][[1]]) else nrow(object$x)
  n.q <- object$n.q
  nu <- object$nu
  if (length(nu) == 1L) {
    nu <- rep(nu, 2L)
  }
  nu.vec <- c(nu[1], rep(nu[2], ncol(object$x.tm) - 1L))
  time.grid <- object$time.unique
  x.tm.grid <- object$x.tm
  l.pred.db.hat.temp <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n), function(i) {
      lapply(seq_len(m.opt[q]), function(m) {
        gamma.step <- gamma[[q]][[m]]
        gamma.i <- boostmtree.predict.match.gamma(gamma.step, membership[[q]][[m]][i])
        as.vector(x.tm.grid %*% c(gamma.i[1, ] * nu.vec))
      })
    })
  })
  l.pred.db.hat <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n), function(i) {
      temp.path <- l.pred.db.hat.temp[[q]][[i]]
      if (identical(object$family, "ordinal") && q > 1L) {
        for (m in seq_len(m.opt[q])) {
          temp.path[[m]] <- pmax(
            temp.path[[m]],
            l.pred.db.hat.temp[[q - 1L]][[i]][[m]]
          )
        }
      }
      Reduce("+", temp.path)
    })
  })
  l.pred.ref.hat <- if (identical(object$family, "nominal")) {
    lapply(seq_len(n), function(i) {
      log((1 + Reduce("+", lapply(seq_len(n.q), function(q) {
        exp(l.pred.db.hat[[q]][[i]])
      })))^(-1))
    })
  } else {
    lapply(seq_len(n), function(i) rep(0, nrow(x.tm.grid)))
  }
  l.pred.hat <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n), function(i) {
      (l.pred.db.hat[[q]][[i]] + l.pred.ref.hat[[i]]) * object$y.sd + object$y.mean
    })
  })
  muhat <- lapply(seq_len(n.q), function(q) {
    lapply(seq_len(n), function(i) {
      boostmtree.get.mu(l.pred.hat[[q]][[i]], family = object$family)
    })
  })
  list(
    time.grid = time.grid,
    muhat = muhat,
    prob.hat.class = boostmtree.build.prob.class(
      mu = muhat,
      family = object$family,
      y.levels = object$y.levels,
      q.set.index = boostmtree.predict.q.set.index(object)
    )
  )
}
generic.predict.boostmtree <- function(
  object,
  x,
  tm,
  id,
  y,
  M = NULL,
  eps = 1e-5,
  use.cv.flag = FALSE,
  partial = FALSE,
  ...
) {
  if (missing(object)) {
    stop("`object` is missing.")
  }
  if (!all(c("boostmtree", "grow") %in% class(object))) {
    stop("`object` must inherit from `c(\"boostmtree\", \"grow\")`.")
  }
  partial <- isTRUE(partial)
  use.cv.flag <- isTRUE(use.cv.flag)
  x.supplied <- !missing(x)
  tm.supplied <- !missing(tm)
  id.supplied <- !missing(id)
  y.supplied <- !missing(y)
  training.names <- object$x.var.names
  na.action <- object$na.action
  if (partial && !x.supplied) {
    stop("`x` must be supplied when `partial = TRUE`.")
  }
  if (partial && !tm.supplied) {
    stop("`tm` must be supplied when `partial = TRUE`.")
  }
  test.flag <- FALSE
  if (!x.supplied) {
    x.subject <- object$x
    time.by.subject <- object$time
    id.long <- object$id
    id.unique <- object$id.unique
    y.raw <- NULL
  } else if (partial) {
    x.subject <- boostmtree.predict.validate.x(x, training.names = training.names)
    time.grid <- sort(unique(tm))
    time.by.subject <- lapply(seq_len(nrow(x.subject)), function(i) time.grid)
    id.unique <- seq_len(nrow(x.subject))
    id.long <- rep(id.unique, each = length(time.grid))
    y.raw <- NULL
  } else if (!tm.supplied || !id.supplied) {
    if (y.supplied) {
      warning("`y` is ignored unless both `tm` and `id` are supplied for the prediction data.")
    }
    x.subject <- boostmtree.predict.validate.x(x, training.names = training.names)
    time.grid <- object$time.unique
    time.by.subject <- lapply(seq_len(nrow(x.subject)), function(i) time.grid)
    id.unique <- seq_len(nrow(x.subject))
    id.long <- rep(id.unique, each = length(time.grid))
    y.raw <- NULL
  } else {
    prepared <- boostmtree.predict.prepare.long.data(
      x = x,
      tm = tm,
      id = id,
      y = if (y.supplied) y else NULL,
      na.action = na.action,
      training.names = training.names
    )
    x.subject <- prepared$x
    time.by.subject <- prepared$time
    id.long <- prepared$id
    id.unique <- prepared$id.unique
    y.raw <- prepared$y
    test.flag <- y.supplied
  }
  n <- nrow(x.subject)
  ni <- vapply(time.by.subject, length, integer(1))
  time.design <- boostmtree.predict.time.design(time.by.subject, object = object)
  prediction.time.unique <- sort(unique(unlist(time.by.subject)))
  if (!is.null(M)) {
    if (length(M) != 1L || is.na(M)) {
      stop("`M` must be NULL or a single positive integer.")
    }
    M <- as.integer(M)
    M <- max(1L, min(M, object$M))
  }
  if (use.cv.flag && !isTRUE(object$cv.flag)) {
    warning("`use.cv.flag = TRUE` requires a fitted object with `cv.flag = TRUE`; using the ordinary fitted coefficients instead.")
    use.cv.flag <- FALSE
  }
  if (use.cv.flag && x.supplied && !partial) {
    warning("`use.cv.flag = TRUE` is only used for the original training subjects; using the ordinary fitted coefficients for new prediction data.")
    use.cv.flag <- FALSE
  }
  if (use.cv.flag && partial && x.supplied && nrow(x.subject) != nrow(object$x)) {
    warning("`use.cv.flag = TRUE` with `partial = TRUE` requires `x` to contain the full fitted subject set in the original row order; using the ordinary fitted coefficients instead.")
    use.cv.flag <- FALSE
  }
  gamma <- if (use.cv.flag) object$gamma.i.list else object$gamma
  membership <- lapply(seq_len(object$n.q), function(q) vector("list", object$M))
  oob.list <- if (use.cv.flag) lapply(seq_len(object$n.q), function(q) vector("list", object$M)) else NULL
  response.info <- NULL
  if (test.flag) {
    response.info <- boostmtree.predict.encode.response(
      y = y.raw,
      object = object,
      id = id.long,
      id.unique = id.unique
    )
  }
  m.max <- if (is.null(M)) object$M else M
  for (q in seq_len(object$n.q)) {
    for (m in seq_len(m.max)) {
      if (!use.cv.flag) {
        prediction.object <- predict.rfsrc(
          object$base.learner[[q]][[m]],
          newdata = x.subject,
          membership = TRUE,
          ptn.count = object$k,
          na.action = object$na.action,
          importance = "none"
        )
        membership[[q]][[m]] <- c(prediction.object$ptn.membership)
      } else {
        oob <- which(object$base.learner[[q]][[m]]$inbag == 0)
        oob.list[[q]][[m]] <- oob
        if (length(oob) > 0L) {
          prediction.object <- predict.rfsrc(
            object$base.learner[[q]][[m]],
            newdata = x.subject[oob, , drop = FALSE],
            membership = TRUE,
            ptn.count = object$k,
            na.action = object$na.action,
            importance = "none"
          )
          membership[[q]][[m]] <- c(prediction.object$ptn.membership)
        } else {
          membership[[q]][[m]] <- character(0)
        }
      }
    }
  }
  l.pred.db.list <- lapply(seq_len(object$n.q), function(q) {
    lapply(seq_len(m.max), function(m) {
      lapply(seq_len(n), function(i) rep(0, ni[i]))
    })
  })
  nu <- object$nu
  if (length(nu) == 1L) {
    nu <- rep(nu, 2L)
  }
  df.time.design <- ncol(object$x.tm)
  nu.vec <- c(nu[1], rep(nu[2], df.time.design - 1L))
  for (q in seq_len(object$n.q)) {
    for (m in seq_len(m.max)) {
      if (!use.cv.flag) {
        gamma.step <- gamma[[q]][[m]]
        gamma.by.subject <- boostmtree.predict.match.gamma(
          gamma.step = gamma.step,
          membership.label = membership[[q]][[m]]
        )
        previous <- if (m == 1L) {
          lapply(seq_len(n), function(i) rep(0, ni[i]))
        } else {
          l.pred.db.list[[q]][[m - 1L]]
        }
        l.pred.db.list[[q]][[m]] <- lapply(seq_len(n), function(i) {
          updated <- previous[[i]] + c(time.design[[i]] %*% c(gamma.by.subject[i, ] * nu.vec))
          if (identical(object$family, "ordinal") && q > 1L) {
            pmax(updated, l.pred.db.list[[q - 1L]][[m]][[i]])
          } else {
            updated
          }
        })
      } else {
        previous <- if (m == 1L) {
          lapply(seq_len(n), function(i) rep(0, ni[i]))
        } else {
          l.pred.db.list[[q]][[m - 1L]]
        }
        l.pred.db.list[[q]][[m]] <- lapply(seq_len(n), function(i) {
          updated <- previous[[i]]
          if (i %in% oob.list[[q]][[m]]) {
            gamma.step <- gamma[[q]][[m]][[i]]
            membership.index <- match(i, oob.list[[q]][[m]])
            gamma.i <- boostmtree.predict.match.gamma(
              gamma.step = gamma.step,
              membership.label = membership[[q]][[m]][membership.index]
            )
            updated <- updated + c(time.design[[i]] %*% c(gamma.i[1, ] * nu.vec))
          }
          if (identical(object$family, "ordinal") && q > 1L) {
            pmax(updated, l.pred.db.list[[q - 1L]][[m]][[i]])
          } else {
            updated
          }
        })
      }
    }
  }
  l.pred.ref <- boostmtree.predict.reference.path(
    l.pred.db.list = l.pred.db.list,
    family = object$family,
    q.total = object$q.total,
    n = n,
    ni = ni,
    m.max = m.max
  )
  l.pred.list <- lapply(seq_len(object$n.q), function(q) {
    lapply(seq_len(m.max), function(m) {
      lapply(seq_len(n), function(i) {
        (l.pred.db.list[[q]][[m]][[i]] + l.pred.ref[[m]][[i]]) * object$y.sd + object$y.mean
      })
    })
  })
  mu.list <- lapply(seq_len(object$n.q), function(q) {
    lapply(seq_len(m.max), function(m) {
      lapply(seq_len(n), function(i) {
        boostmtree.get.mu(l.pred.list[[q]][[m]][[i]], family = object$family)
      })
    })
  })
  err.rate.raw <- NULL
  if (test.flag) {
    err.rate.raw <- lapply(seq_len(object$n.q), function(q) {
      err <- matrix(
        vapply(seq_len(m.max), function(m) {
          c(
            boostmtree.l1.dist(response.info$y.org[[q]], mu.list[[q]][[m]]),
            boostmtree.l2.dist(response.info$y.org[[q]], mu.list[[q]][[m]])
          )
        }, numeric(2)),
        ncol = 2L,
        byrow = TRUE
      )
      colnames(err) <- c("l1", "l2")
      err
    })
  }
  m.opt <- boostmtree.predict.select.m(
    object = object,
    err.rate = err.rate.raw,
    M = M,
    eps = eps
  )
  mu <- lapply(seq_len(object$n.q), function(q) {
    mu.list[[q]][[m.opt[q]]]
  })
  prob.class <- boostmtree.build.prob.class(
    mu = mu,
    family = object$family,
    y.levels = object$y.levels,
    q.set.index = boostmtree.predict.q.set.index(object)
  )
  extrapolated <- if (!use.cv.flag) {
    boostmtree.predict.extrapolate(
      object = object,
      membership = membership,
      gamma = object$gamma,
      m.opt = m.opt
    )
  } else {
    list(
      time.grid = NULL,
      muhat = NULL,
      prob.hat.class = NULL
    )
  }
  rmse <- if (test.flag) {
    vapply(seq_len(object$n.q), function(q) err.rate.raw[[q]][m.opt[q], "l2"], numeric(1))
  } else {
    NULL
  }
  err.rate <- if (!is.null(err.rate.raw)) {
    lapply(seq_len(object$n.q), function(q) err.rate.raw[[q]] / object$y.sd)
  } else {
    NULL
  }
  boost.obj <- object
  boost.obj$base.learner <- NULL
  boost.obj$membership <- NULL
  boost.obj$gamma <- NULL
  boost.obj$gamma.i.list <- NULL
  predict.object <- list(
    boost.obj = boost.obj,
    x = x.subject,
    time = time.by.subject,
    time.unique = prediction.time.unique,
    time.grid = extrapolated$time.grid,
    id = id.long,
    id.unique = id.unique,
    y = if (test.flag) response.info$y else NULL,
    y.org = if (test.flag) {
      boostmtree.predict.flatten.single.response(response.info$y.org, family = object$family)
    } else {
      NULL
    },
    family = object$family,
    y.mean = object$y.mean,
    y.sd = object$y.sd,
    y.levels = object$y.levels,
    y.reference = object$y.reference,
    x.var.names = object$x.var.names,
    k = object$k,
    n = n,
    ni = ni,
    n.q = object$n.q,
    q.total = object$q.total,
    q.set = object$q.set,
    nu = object$nu,
    nu.vec = nu.vec,
    time.design = time.design,
    df.time.design = ncol(object$x.tm),
    base.learner = object$base.learner,
    gamma = if (use.cv.flag) object$gamma.i.list else object$gamma,
    membership = membership,
    mu = boostmtree.predict.flatten.single.response(mu, family = object$family),
    prob.class = prob.class,
    muhat = boostmtree.predict.flatten.single.response(extrapolated$muhat, family = object$family),
    prob.hat.class = extrapolated$prob.hat.class,
    err.rate = boostmtree.predict.flatten.error(err.rate, family = object$family),
    rmse = if (!is.null(rmse)) as.numeric(rmse / object$y.sd) else NULL,
    m.opt = m.opt,
    partial = partial,
    use.cv.flag = use.cv.flag
  )
  class(predict.object) <- c("boostmtree", "predict", class(object)[3])
  invisible(predict.object)
}
