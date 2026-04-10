`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
boostmtree.is.grow.object <- function(x) {
  inherits(x, "boostmtree") && identical(class(x)[2], "grow")
}
boostmtree.is.predict.object <- function(x) {
  inherits(x, "boostmtree") && (
    identical(class(x)[2], "predict") || !is.null(x$boost.obj)
  )
}
boostmtree.as.q.list <- function(x, n.q) {
  if (is.null(x)) {
    return(NULL)
  }
  if (n.q == 1L) {
    if (!is.list(x)) {
      return(list(x))
    }
    if (length(x) > 0L && !is.list(x[[1]])) {
      return(list(x))
    }
  }
  x
}
boostmtree.as.error.list <- function(x, n.q) {
  if (is.null(x)) {
    return(NULL)
  }
  if (n.q == 1L && is.matrix(x)) {
    return(list(x))
  }
  x
}
boostmtree.extract.iteration.value <- function(path, m.opt, q) {
  if (is.null(path) || is.null(m.opt) || length(m.opt) < q) {
    return(NULL)
  }
  iteration.index <- as.integer(m.opt[q])
  if (length(iteration.index) != 1L || is.na(iteration.index) || iteration.index < 1L) {
    return(NULL)
  }
  if (is.matrix(path)) {
    if (nrow(path) < iteration.index || ncol(path) < q) {
      return(NULL)
    }
    return(path[iteration.index, q])
  }
  if (is.atomic(path) && !is.list(path)) {
    if (length(path) < iteration.index) {
      return(NULL)
    }
    return(path[iteration.index])
  }
  NULL
}
boostmtree.count.unique.time <- function(time.by.subject) {
  if (is.null(time.by.subject)) {
    return(0L)
  }
  length(sort(unique(unlist(time.by.subject))))
}
boostmtree.resolve.display.object <- function(x) {
  ## Predict objects must be detected first and explicitly.
  ## They contain a fitted grow object in `boost.obj`, but the display
  ## should use prediction-side summaries and labels.
  if (boostmtree.is.predict.object(x)) {
    if (is.null(x$boost.obj)) {
      stop("Predict objects must contain a `boost.obj` component.")
    }
    model <- x$boost.obj
    return(list(
      kind = "predict",
      model = model,
      fit = x,
      x = x$x %||% model$x,
      time = x$time %||% model$time,
      mu = boostmtree.as.q.list(x$mu %||% model$mu, model$n.q),
      y.org = boostmtree.as.q.list(x$y.org, model$n.q),
      err.rate = boostmtree.as.error.list(x$err.rate, model$n.q),
      rmse = x$rmse,
      m.opt = x$m.opt %||% model$m.opt,
      vimp = x$vimp %||% NULL
    ))
  }
  if (boostmtree.is.grow.object(x)) {
    model <- x
    return(list(
      kind = "grow",
      model = model,
      fit = x,
      x = x$x,
      time = x$time,
      mu = boostmtree.as.q.list(x$mu, model$n.q),
      y.org = boostmtree.as.q.list(x$y.org, model$n.q),
      err.rate = boostmtree.as.error.list(x$err.rate, model$n.q),
      rmse = x$rmse,
      m.opt = x$m.opt,
      vimp = NULL
    ))
  }
  stop("This method only works for `(boostmtree, grow)` or `(boostmtree, predict)` objects.")
}
boostmtree.plot.data <- function(x, use.rmse = TRUE) {
  display <- boostmtree.resolve.display.object(x)
  model <- display$model
  plot.data <- lapply(seq_len(model$n.q), function(q) {
    error.path <- NULL
    if (!is.null(display$err.rate)) {
      error.value <- display$err.rate[[q]][, "l2"]
      if (!use.rmse) {
        error.value <- (error.value * model$y.sd)^2
      }
      error.path <- data.frame(
        iteration = seq_len(if (is.matrix(display$err.rate[[q]])) nrow(display$err.rate[[q]]) else length(error.value)),
        value = error.value
      )
    }
    list(
      q.index = q,
      q.label = if (!is.null(model$q.set) && !all(is.na(model$q.set))) model$q.set[q] else NULL,
      time.by.subject = display$time,
      fitted.by.subject = display$mu[[q]],
      observed.by.subject = if (!is.null(display$y.org)) display$y.org[[q]] else NULL,
      error.path = error.path,
      rho.path = if (!is.null(model$rho)) {
        if (is.matrix(model$rho)) data.frame(iteration = seq_len(nrow(model$rho)), value = model$rho[, q]) else data.frame(iteration = seq_along(model$rho), value = model$rho)
      } else {
        NULL
      },
      phi.path = if (!is.null(model$phi)) {
        if (is.matrix(model$phi)) data.frame(iteration = seq_len(nrow(model$phi)), value = model$phi[, q]) else data.frame(iteration = seq_along(model$phi), value = model$phi)
      } else {
        NULL
      },
      lambda.path = if (!is.null(model$lambda)) {
        if (is.matrix(model$lambda)) data.frame(iteration = seq_len(nrow(model$lambda)), value = model$lambda[, q]) else data.frame(iteration = seq_along(model$lambda), value = model$lambda)
      } else {
        NULL
      },
      variable.importance = display$vimp,
      m.opt = display$m.opt %||% NULL
    )
  })
  structure(
    plot.data,
    family = model$family,
    kind = display$kind
  )
}
# Backward-compatible internal alias for plot helper during refactor.
boostmtree.plot.payload <- function(x, use.rmse = TRUE) {
  boostmtree.plot.data(x, use.rmse = use.rmse)
}
