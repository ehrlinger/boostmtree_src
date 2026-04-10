marginal.plot <- function(object, ...) {
  UseMethod("marginal.plot")
}
marginal.plot.boostmtree <- function(
  object,
  M = NULL,
  x.var.names = NULL,
  time.points = NULL,
  subset = NULL,
  prob.class = FALSE,
  response.labels = NULL,
  output = c("plot", "data", "pdf"),
  file = NULL,
  verbose = TRUE,
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
  if (is.null(M)) {
    M <- if (!is.null(object$m.opt)) max(object$m.opt, na.rm = TRUE) else object$M
  } else {
    M <- max(1L, min(as.integer(M), object$M))
  }
  pred.object <- predict(object, M = M, ...)
  pred.response <- lapply(seq_len(response.info$n.response), function(q.out) {
    q.index <- response.info$response.index[q.out]
    boostmtree.effect.extract.marginal.response(
      pred.object = pred.object,
      object = object,
      q.index = q.index,
      prob.class = response.info$prob.class
    )
  })
  if (!is.null(subset)) {
    pred.response <- lapply(pred.response, function(q.list) q.list[subset])
  }
  work.x <- boostmtree.effect.subset.data(object$x, subset = subset)
  raw.data <- smooth.data <- vector("list", response.info$n.response)
  names(raw.data) <- names(smooth.data) <- response.info$response.labels
  for (q in seq_len(response.info$n.response)) {
    response.raw <- response.smooth <- vector("list", length(x.var.names))
    names(response.raw) <- names(response.smooth) <- x.var.names
    for (nm in x.var.names) {
      x.column <- work.x[[nm]]
      raw.by.time <- smooth.by.time <- vector("list", length(time.info$index))
      names(raw.by.time) <- names(smooth.by.time) <- paste("time =", format(signif(time.info$time.points, 4), trim = TRUE))
      for (j in seq_along(time.info$index)) {
        y.value <- vapply(seq_along(pred.response[[q]]), function(i) {
          pred.response[[q]][[i]][time.info$index[j]]
        }, numeric(1))
        raw.by.time[[j]] <- data.frame(
          x = x.column,
          y = y.value,
          stringsAsFactors = FALSE
        )
        smooth.by.time[[j]] <- boostmtree.effect.smooth.xy(x.column, y.value)
      }
      response.raw[[nm]] <- raw.by.time
      response.smooth[[nm]] <- smooth.by.time
    }
    raw.data[[q]] <- response.raw
    smooth.data[[q]] <- response.smooth
  }
  result <- list(
    data = boostmtree.effect.flatten.single.response(raw.data, response.info$n.response),
    smooth = boostmtree.effect.flatten.single.response(smooth.data, response.info$n.response),
    time.points = time.info$time.points,
    x.var.names = x.var.names,
    response.labels = response.info$response.labels,
    family = object$family,
    prob.class = response.info$prob.class,
    M = M,
    call = match.call()
  )
  class(result) <- c("marginal.plot.boostmtree", "boostmtree.effect.plot")
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
plot.marginal.plot.boostmtree <- function(x, output = c("plot", "pdf"), file = NULL, verbose = TRUE, ...) {
  output <- match.arg(output)
  if (identical(output, "pdf")) {
    file <- boostmtree.effect.device.file(file, "marginal_plot_boostmtree.pdf")
    grDevices::pdf(file = file, width = 10, height = 10)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  smooth.list <- if (length(x$response.labels) == 1L) {
    stats::setNames(list(x$smooth), x$response.labels)
  } else {
    x$smooth
  }
  response.names <- if (is.null(names(smooth.list))) x$response.labels else names(smooth.list)
  n.panel <- sum(vapply(smooth.list, length, integer(1)))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par), add = TRUE)
  graphics::par(mfrow = boostmtree.effect.layout(n.panel))
  for (response.index in seq_along(response.names)) {
    panel.list <- smooth.list[[response.index]]
    for (nm in names(panel.list)) {
      boostmtree.effect.draw.curves(
        curve.list = panel.list[[nm]],
        time.points = x$time.points,
        x.label = nm,
        y.label = if (isTRUE(x$prob.class)) "predicted probability" else "predicted response",
        main = if (length(response.names) > 1L) paste("response:", response.names[response.index]) else NULL
      )
    }
  }
  if (identical(output, "pdf") && isTRUE(verbose)) {
    cat("Plot saved at:", file, "\n")
  }
  invisible(x)
}
