boostmtree.draw.single.panel <- function(panel, model, payload.q, use.rmse) {
  if (identical(panel, "trajectory")) {
    plot(
      unlist(payload.q$time.by.subject),
      unlist(payload.q$fitted.by.subject),
      xlab = "time",
      ylab = "fitted",
      type = "n"
    )
    boostmtree.line.plot(payload.q$time.by.subject, payload.q$fitted.by.subject)
    return(invisible(NULL))
  }
  if (identical(panel, "residual")) {
    residual.by.subject <- lapply(seq_along(payload.q$observed.by.subject), function(i) {
      payload.q$observed.by.subject[[i]] - payload.q$fitted.by.subject[[i]]
    })
    plot(
      unlist(payload.q$time.by.subject),
      unlist(residual.by.subject),
      xlab = "time",
      ylab = "residual",
      type = "n"
    )
    boostmtree.line.plot(payload.q$time.by.subject, residual.by.subject)
    abline(h = 0, lty = 2, col = "gray")
    return(invisible(NULL))
  }
  if (identical(panel, "observed.vs.fitted")) {
    plot(
      unlist(payload.q$observed.by.subject),
      unlist(payload.q$fitted.by.subject),
      xlab = "observed",
      ylab = "fitted",
      type = "n"
    )
    if (isTRUE(model$univariate)) {
      boostmtree.point.plot(payload.q$observed.by.subject, payload.q$fitted.by.subject)
    } else {
      boostmtree.line.plot(payload.q$observed.by.subject, payload.q$fitted.by.subject)
    }
    abline(0, 1, col = "gray", lty = 2)
    return(invisible(NULL))
  }
  if (identical(panel, "error.path")) {
    y.label <- if (use.rmse) "standardized RMSE" else "MSE"
    plot(
      payload.q$error.path$iteration,
      payload.q$error.path$value,
      xlab = "iteration",
      ylab = y.label,
      type = "l",
      lty = 1
    )
    if (!is.null(payload.q$m.opt)) {
      abline(v = payload.q$m.opt[payload.q$q.index], lty = 2, col = 2, lwd = 2)
    }
    return(invisible(NULL))
  }
  if (identical(panel, "rho.path")) {
    smooth.path <- boostmtree.lowess(payload.q$rho.path$iteration, payload.q$rho.path$value, f = 0.5)
    plot(
      payload.q$rho.path$iteration,
      payload.q$rho.path$value,
      xlab = "iteration",
      ylab = expression(rho),
      type = "n"
    )
    lines(smooth.path)
    if (!is.null(payload.q$m.opt)) {
      abline(v = payload.q$m.opt[payload.q$q.index], lty = 2, col = 2, lwd = 2)
    }
    return(invisible(NULL))
  }
  if (identical(panel, "phi.path")) {
    smooth.path <- boostmtree.lowess(payload.q$phi.path$iteration, payload.q$phi.path$value, f = 0.5)
    plot(
      payload.q$phi.path$iteration,
      payload.q$phi.path$value,
      xlab = "iteration",
      ylab = expression(phi),
      type = "n"
    )
    lines(smooth.path)
    if (!is.null(payload.q$m.opt)) {
      abline(v = payload.q$m.opt[payload.q$q.index], lty = 2, col = 2, lwd = 2)
    }
    return(invisible(NULL))
  }
  if (identical(panel, "lambda.path")) {
    smooth.path <- boostmtree.lowess(payload.q$lambda.path$iteration, payload.q$lambda.path$value, f = 0.5)
    plot(
      payload.q$lambda.path$iteration,
      payload.q$lambda.path$value,
      xlab = "iteration",
      ylab = expression(lambda),
      type = "n"
    )
    lines(smooth.path)
    if (!is.null(payload.q$m.opt)) {
      abline(v = payload.q$m.opt[payload.q$q.index], lty = 2, col = 2, lwd = 2)
    }
    return(invisible(NULL))
  }
  if (identical(panel, "variable.importance")) {
    barplot(payload.q$variable.importance, las = 2, ylab = "variable importance")
    return(invisible(NULL))
  }
  plot.new()
  title(panel)
  invisible(NULL)
}
boostmtree.draw.single.plot <- function(payload.q, model, use.rmse = TRUE) {
  panels <- character(0)
  if (!is.null(payload.q$time.by.subject) && !is.null(payload.q$fitted.by.subject) && !isTRUE(model$univariate)) {
    panels <- c(panels, "trajectory")
  }
  if (!is.null(payload.q$observed.by.subject) && !is.null(payload.q$fitted.by.subject)) {
    if (!is.null(payload.q$error.path)) {
      panels <- c(panels, "observed.vs.fitted")
    } else {
      panels <- c(panels, "residual", "observed.vs.fitted")
    }
  }
  if (!is.null(payload.q$error.path)) {
    panels <- c(panels, "error.path")
  }
  if (!is.null(payload.q$rho.path)) {
    panels <- c(panels, "rho.path")
  }
  if (!is.null(payload.q$phi.path)) {
    panels <- c(panels, "phi.path")
  }
  if (!is.null(payload.q$lambda.path)) {
    panels <- c(panels, "lambda.path")
  }
  if (!is.null(payload.q$variable.importance)) {
    panels <- c(panels, "variable.importance")
  }
  if (length(panels) == 0L) {
    panels <- "trajectory"
  }
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par), add = TRUE)
  n.panel <- length(panels)
  if (n.panel == 1L) {
    par(mfrow = c(1, 1))
  } else if (n.panel == 2L) {
    par(mfrow = c(1, 2))
  } else {
    par(mfrow = c(ceiling(n.panel / 2), 2))
  }
  for (panel in panels) {
    boostmtree.draw.single.panel(
      panel = panel,
      model = model,
      payload.q = payload.q,
      use.rmse = use.rmse
    )
  }
  if (!is.null(payload.q$q.label)) {
    mtext(
      side = 3,
      line = -1.5,
      outer = FALSE,
      text = paste("response level:", payload.q$q.label)
    )
  }
  invisible(NULL)
}
plot.boostmtree <- function(
  x,
  use.rmse = TRUE,
  output = c("plot", "data", "pdf"),
  file = NULL,
  width = 10,
  height = 10,
  verbose = TRUE,
  ...
) {
  output <- match.arg(output)
  display <- boostmtree.resolve.display.object(x)
  payload <- boostmtree.plot.payload(x, use.rmse = use.rmse)
  if (identical(output, "data")) {
    return(payload)
  }
  if (identical(output, "pdf")) {
    if (is.null(file)) {
      stop("`file` must be supplied when `output = \"pdf\"`.")
    }
    pdf(file = file, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
  }
  for (q in seq_along(payload)) {
    boostmtree.draw.single.plot(
      payload.q = payload[[q]],
      model = display$model,
      use.rmse = use.rmse
    )
  }
  if (identical(output, "pdf") && isTRUE(verbose)) {
    cat("Plot saved to:", file, "\n")
  }
  invisible(payload)
}
