# Shared helpers for effect-plot summaries ---------------------------------
boostmtree.effect.is.grow.object <- function(object) {
  inherits(object, c("boostmtree", "grow"))
}
boostmtree.effect.as.q.list <- function(x, n.q) {
  if (is.null(x)) {
    return(NULL)
  }
  if (n.q == 1L) {
    if (!is.list(x)) {
      return(list(x))
    }
    if (length(x) == 0L) {
      return(list(x))
    }
    first <- x[[1L]]
    if (!is.list(first)) {
      return(list(x))
    }
  }
  x
}
boostmtree.effect.flatten.single.response <- function(x, n.q) {
  if (is.null(x)) {
    return(NULL)
  }
  if (n.q == 1L) {
    return(x[[1L]])
  }
  x
}
boostmtree.effect.match.x.var.names <- function(object, x.var.names) {
  if (is.null(x.var.names)) {
    x.var.names <- object$x.var.names
  }
  x.var.names <- intersect(x.var.names, object$x.var.names)
  if (length(x.var.names) == 0L) {
    stop("`x.var.names` does not match any fitted covariate names.")
  }
  x.var.names
}
boostmtree.effect.time.info <- function(object, time.points = NULL) {
  observed.time <- if (!is.null(object$time.unique)) {
    object$time.unique
  } else {
    sort(unique(unlist(object$time)))
  }
  if (is.null(time.points)) {
    chosen <- unique(quantile(observed.time, probs = (1:9) / 10, na.rm = TRUE))
  } else {
    chosen <- time.points
  }
  index <- vapply(chosen, function(tt) {
    which.min(abs(observed.time - tt))
  }, integer(1))
  actual <- observed.time[index]
  column.names <- make.names(paste0("time.", format(signif(actual, 4), trim = TRUE)))
  list(
    observed.time = observed.time,
    index = index,
    time.points = actual,
    column.names = column.names
  )
}
boostmtree.effect.response.info <- function(object, prob.class = FALSE) {
  family <- object$family
  if (prob.class && !identical(family, "ordinal")) {
    warning("`prob.class = TRUE` is only used for the ordinal family; using fitted means instead.")
    prob.class <- FALSE
  }
  if (prob.class) {
    response.labels <- as.character(object$y.levels)
    n.response <- length(response.labels)
  } else if (identical(family, "continuous")) {
    response.labels <- "response"
    n.response <- 1L
  } else if (identical(family, "binary")) {
    response.labels <- as.character(object$q.set[1])
    n.response <- 1L
  } else {
    response.labels <- as.character(object$q.set)
    n.response <- length(response.labels)
  }
  list(
    prob.class = prob.class,
    response.labels = response.labels,
    n.response = n.response
  )
}
boostmtree.effect.select.responses <- function(response.info, response.labels = NULL) {
  response.info$response.index <- seq_len(response.info$n.response)
  if (is.null(response.labels)) {
    return(response.info)
  }
  response.labels <- as.character(response.labels)
  matched <- match(response.labels, response.info$response.labels)
  if (anyNA(matched)) {
    stop(
      "`response.labels` must match one or more available response labels: ",
      paste(response.info$response.labels, collapse = ", ")
    )
  }
  matched <- unique(matched)
  response.info$response.index <- matched
  response.info$response.labels <- response.info$response.labels[matched]
  response.info$n.response <- length(matched)
  response.info
}
boostmtree.effect.subset.data <- function(x, subset) {
  if (is.null(subset)) {
    return(x)
  }
  x[subset, , drop = FALSE]
}
boostmtree.effect.validate.conditioning <- function(object, x.data, conditional.x.var.names, conditional.values) {
  if (is.null(conditional.x.var.names) && is.null(conditional.values)) {
    return(x.data)
  }
  if (is.null(conditional.x.var.names) || is.null(conditional.values)) {
    stop("`conditional.x.var.names` and `conditional.values` must be supplied together.")
  }
  if (length(conditional.x.var.names) != length(conditional.values)) {
    stop("`conditional.x.var.names` and `conditional.values` must have the same length.")
  }
  missing.names <- setdiff(conditional.x.var.names, colnames(x.data))
  if (length(missing.names) > 0L) {
    stop(
      "Conditioning variables not found in the fitted data: ",
      paste(missing.names, collapse = ", ")
    )
  }
  out <- x.data
  for (i in seq_along(conditional.x.var.names)) {
    nm <- conditional.x.var.names[i]
    value <- conditional.values[i]
    template <- out[[nm]]
    if (is.factor(template)) {
      level.values <- levels(template)
      if (!as.character(value) %in% level.values) {
        stop(
          "Conditioning value `", value, "` is not an observed level of `", nm, "`."
        )
      }
      out[[nm]] <- factor(rep(as.character(value), nrow(out)), levels = level.values)
    } else {
      out[[nm]] <- rep(value, nrow(out))
    }
  }
  out
}
boostmtree.effect.resolve.x.values <- function(x.column, x.var.values, x.var.name, n.points) {
  if (!is.null(x.var.values)) {
    if (!is.list(x.var.values)) {
      stop("`x.var.values` must be NULL or a named list.")
    }
    if (is.null(names(x.var.values)) || !x.var.name %in% names(x.var.values)) {
      stop("`x.var.values` must be a named list containing `", x.var.name, "`.")
    }
    values <- x.var.values[[x.var.name]]
  } else {
    if (is.factor(x.column)) {
      values <- levels(droplevels(x.column))
    } else {
      observed <- sort(unique(stats::na.omit(x.column)))
      if (length(observed) == 0L) {
        stop("No observed values are available for `", x.var.name, "`.")
      }
      pick <- unique(as.integer(seq(1, length(observed), length.out = min(n.points, length(observed)))))
      values <- observed[pick]
    }
  }
  if (length(values) == 0L) {
    stop("No values are available for `", x.var.name, "`.")
  }
  list(
    values = values,
    is.factor = is.factor(x.column),
    display = if (is.factor(x.column)) as.character(values) else values
  )
}
boostmtree.effect.assign.column <- function(template, values) {
  if (is.factor(template)) {
    factor(rep(as.character(values), length.out = length(template)), levels = levels(template))
  } else if (is.character(template)) {
    rep(as.character(values), length.out = length(template))
  } else {
    rep(values, length.out = length(template))
  }
}
boostmtree.effect.extract.partial.response <- function(pred.object, object, q.index, prob.class = FALSE) {
  if (prob.class) {
    return(pred.object$prob.class[[q.index]])
  }
  boostmtree.effect.as.q.list(pred.object$mu, if (identical(object$family, "continuous") || identical(object$family, "binary")) 1L else object$n.q)[[q.index]]
}
boostmtree.effect.extract.marginal.response <- function(pred.object, object, q.index, prob.class = FALSE) {
  if (prob.class) {
    return(pred.object$prob.hat.class[[q.index]])
  }
  boostmtree.effect.as.q.list(pred.object$muhat, if (identical(object$family, "continuous") || identical(object$family, "binary")) 1L else object$n.q)[[q.index]]
}
boostmtree.effect.normalize.response.list <- function(x, response.labels) {
  if (length(response.labels) == 1L) {
    return(list(setNames(list(x), response.labels))[[1L]])  # not used directly
  }
  x
}
boostmtree.effect.curves.as.response.list <- function(x, response.labels) {
  if (length(response.labels) == 1L) {
    list(setNames(list(x), response.labels))[[1L]]
  } else {
    x
  }
}
boostmtree.effect.make.plot.data <- function(curve.data.frame) {
  list(
    x = curve.data.frame[[1L]],
    y.list = lapply(seq_len(ncol(curve.data.frame) - 1L), function(j) {
      data.frame(
        x = curve.data.frame[[1L]],
        y = curve.data.frame[[j + 1L]],
        stringsAsFactors = FALSE
      )
    })
  )
}
boostmtree.effect.smooth.xy <- function(x, y) {
  if (is.factor(x) || is.character(x)) {
    x.factor <- factor(x)
    y.mean <- tapply(y, x.factor, mean, na.rm = TRUE)
    data.frame(
      x = names(y.mean),
      y = as.numeric(y.mean),
      stringsAsFactors = FALSE
    )
  } else {
    fit <- if (all(is.na(y)) || stats::sd(y, na.rm = TRUE) == 0) {
      list(x = x, y = y)
    } else {
      stats::lowess(x, y)
    }
    data.frame(x = fit$x, y = fit$y, stringsAsFactors = FALSE)
  }
}
boostmtree.effect.draw.curves <- function(curve.list, time.points, x.label, y.label, main = NULL) {
  first.x <- curve.list[[1L]]$x
  is.discrete <- is.factor(first.x) || is.character(first.x)
  if (is.discrete) {
    x.labels <- as.character(first.x)
    x.pos <- seq_along(x.labels)
    y.range <- range(unlist(lapply(curve.list, function(curve) curve$y)), na.rm = TRUE)
    plot(
      x = x.pos,
      y = curve.list[[1L]]$y,
      type = "n",
      xaxt = "n",
      xlab = x.label,
      ylab = y.label,
      ylim = y.range,
      main = main
    )
    axis(1, at = x.pos, labels = x.labels)
    for (j in seq_along(curve.list)) {
      lines(x.pos, curve.list[[j]]$y, type = "b", col = j, pch = 16)
    }
  } else {
    x.range <- range(unlist(lapply(curve.list, function(curve) curve$x)), na.rm = TRUE)
    y.range <- range(unlist(lapply(curve.list, function(curve) curve$y)), na.rm = TRUE)
    plot(
      x = curve.list[[1L]]$x,
      y = curve.list[[1L]]$y,
      type = "n",
      xlab = x.label,
      ylab = y.label,
      xlim = x.range,
      ylim = y.range,
      main = main
    )
    for (j in seq_along(curve.list)) {
      lines(curve.list[[j]]$x, curve.list[[j]]$y, col = j, lty = 1)
    }
  }
  if (length(curve.list) > 1L) {
    legend(
      "topright",
      legend = format(signif(time.points, 4), trim = TRUE),
      col = seq_along(curve.list),
      lty = 1,
      pch = if (is.discrete) 16 else NA,
      bty = "n",
      title = "time"
    )
  }
  invisible(NULL)
}
boostmtree.effect.device.file <- function(file, default.name) {
  if (is.null(file)) {
    return(file.path(tempdir(), default.name))
  }
  file
}
boostmtree.effect.layout <- function(n.panel) {
  nr <- ceiling(sqrt(n.panel))
  nc <- ceiling(n.panel / nr)
  c(nr, nc)
}
# Partial-plot generic and methods ------------------------------------------
partial.plot <- function(object, ...) {
  UseMethod("partial.plot")
}
partial.plot.boostmtree <- function(
  object,
  M = NULL,
  x.var.names = NULL,
  time.points = NULL,
  x.var.values = NULL,
  n.points = 25,
  subset = NULL,
  prob.class = FALSE,
  response.labels = NULL,
  conditional.x.var.names = NULL,
  conditional.values = NULL,
  output = c("plot", "data", "pdf"),
  file = NULL,
  verbose = TRUE,
  use.cv.flag = NULL,
  ...
) {
  if (!boostmtree.effect.is.grow.object(object)) {
    stop("This function only works for objects of class `(boostmtree, grow)`.")
  }
  output <- match.arg(output)
  x.var.names <- boostmtree.effect.match.x.var.names(object, x.var.names)
  time.info <- boostmtree.effect.time.info(object, time.points)
  response.info <- boostmtree.effect.response.info(object, prob.class = prob.class)
  response.info <- boostmtree.effect.select.responses(response.info, response.labels = response.labels)
  work.x <- boostmtree.effect.subset.data(object$x, subset = subset)
  work.x <- boostmtree.effect.validate.conditioning(
    object = object,
    x.data = work.x,
    conditional.x.var.names = conditional.x.var.names,
    conditional.values = conditional.values
  )
  if (is.null(M)) {
    M <- if (!is.null(object$m.opt)) max(object$m.opt, na.rm = TRUE) else object$M
  } else {
    M <- max(1L, min(as.integer(M), object$M))
  }
  if (is.null(use.cv.flag)) {
    use.cv.flag <- isTRUE(object$cv.flag)
  } else {
    use.cv.flag <- isTRUE(use.cv.flag)
  }
  if (use.cv.flag && !isTRUE(object$cv.flag)) {
    warning("`use.cv.flag = TRUE` requires a fitted object with `cv.flag = TRUE`; using the ordinary fitted coefficients instead.")
    use.cv.flag <- FALSE
  }
  if (use.cv.flag && !is.null(subset)) {
    warning("`use.cv.flag = TRUE` is only supported when all fitted subjects are used; using the ordinary fitted coefficients instead.")
    use.cv.flag <- FALSE
  }
  curves <- vector("list", response.info$n.response)
  names(curves) <- response.info$response.labels
  for (q.out in seq_len(response.info$n.response)) {
    q.index <- response.info$response.index[q.out]
    response.curves <- lapply(x.var.names, function(nm) {
      x.info <- boostmtree.effect.resolve.x.values(
        x.column = work.x[[nm]],
        x.var.values = x.var.values,
        x.var.name = nm,
        n.points = n.points
      )
      curve.matrix <- matrix(
        NA_real_,
        nrow = length(x.info$values),
        ncol = length(time.info$index)
      )
      for (k in seq_along(x.info$values)) {
        new.x <- work.x
        new.x[[nm]] <- boostmtree.effect.assign.column(new.x[[nm]], x.info$values[[k]])
        pred.object <- predict(
          object,
          x = new.x,
          tm = time.info$observed.time,
          M = M,
          partial = TRUE,
          use.cv.flag = use.cv.flag,
          ...
        )
        pred.response <- boostmtree.effect.extract.partial.response(
          pred.object = pred.object,
          object = object,
          q.index = q.index,
          prob.class = response.info$prob.class
        )
        curve.matrix[k, ] <- colMeans(
          do.call(
            rbind,
            lapply(pred.response, function(subject.curve) subject.curve[time.info$index])
          ),
          na.rm = TRUE
        )
      }
      out <- data.frame(
        x = x.info$display,
        curve.matrix,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      names(out) <- c("x", time.info$column.names)
      out
    })
    names(response.curves) <- x.var.names
    curves[[q.out]] <- response.curves
  }
  result <- list(
    curves = boostmtree.effect.flatten.single.response(curves, response.info$n.response),
    time.points = time.info$time.points,
    x.var.names = x.var.names,
    response.labels = response.info$response.labels,
    family = object$family,
    prob.class = response.info$prob.class,
    M = M,
    use.cv.flag = use.cv.flag,
    call = match.call()
  )
  class(result) <- c("partial.plot.boostmtree", "boostmtree.effect.plot")
  if (identical(output, "data")) {
    return(result)
  }
  if (identical(output, "pdf")) {
    plot(result, output = "pdf", file = file, verbose = verbose, ...)
    return(invisible(result))
  }
  plot(result, ...)
  invisible(result)
}
plot.partial.plot.boostmtree <- function(x, output = c("plot", "pdf"), file = NULL, verbose = TRUE, ...) {
  output <- match.arg(output)
  if (identical(output, "pdf")) {
    file <- boostmtree.effect.device.file(file, "partial_plot_boostmtree.pdf")
    grDevices::pdf(file = file, width = 10, height = 10)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  curve.list <- if (length(x$response.labels) == 1L) {
    stats::setNames(list(x$curves), x$response.labels)
  } else {
    x$curves
  }
  response.names <- if (is.null(names(curve.list))) x$response.labels else names(curve.list)
  n.panel <- sum(vapply(curve.list, length, integer(1)))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par), add = TRUE)
  graphics::par(mfrow = boostmtree.effect.layout(n.panel))
  for (response.index in seq_along(response.names)) {
    panel.list <- curve.list[[response.index]]
    for (nm in names(panel.list)) {
      panel.data <- panel.list[[nm]]
      curve.data <- lapply(seq_len(ncol(panel.data) - 1L), function(j) {
        data.frame(
          x = panel.data[[1L]],
          y = panel.data[[j + 1L]],
          stringsAsFactors = FALSE
        )
      })
      boostmtree.effect.draw.curves(
        curve.list = curve.data,
        time.points = x$time.points,
        x.label = nm,
        y.label = if (isTRUE(x$prob.class)) "predicted probability" else "predicted response (adjusted)",
        main = if (length(response.names) > 1L) paste("response:", response.names[response.index]) else NULL
      )
    }
  }
  if (identical(output, "pdf") && isTRUE(verbose)) {
    cat("Plot saved at:", file, "\n")
  }
  invisible(x)
}
