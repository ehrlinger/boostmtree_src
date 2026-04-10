boostmtree.vimp.plot.data <- function(x,
                                      show.interaction = TRUE,
                                      show.time.effect = TRUE,
                                      use.percent = TRUE) {
  if (!inherits(x, "vimp.boostmtree")) {
    stop("`x` must inherit from `vimp.boostmtree`.")
  }
  scale.factor <- if (isTRUE(use.percent)) 100 else 1
  q.names <- colnames(x$main)
  plot.data <- lapply(seq_len(ncol(x$main)), function(q) {
    list(
      q.index = q,
      q.name = q.names[q],
      main = x$main[, q] * scale.factor,
      interaction = if (!is.null(x$interaction) && isTRUE(show.interaction)) x$interaction[, q] * scale.factor else NULL,
      time.effect = if (!is.null(x$time.effect) && isTRUE(show.time.effect)) x$time.effect[q] * scale.factor else NULL,
      x.var.names = rownames(x$main)
    )
  })
  structure(plot.data, use.percent = use.percent)
}
plot.vimp.boostmtree <- function(x,
                                 show.interaction = TRUE,
                                 show.time.effect = TRUE,
                                 output = c("plot", "data", "pdf"),
                                 file = NULL,
                                 main = "Variable importance (%)",
                                 col = grey(0.80),
                                 cex.names = 0.8,
                                 eps = 0.1,
                                 ...) {
  if (!inherits(x, "vimp.boostmtree")) {
    stop("`x` must inherit from `vimp.boostmtree`.")
  }
  output <- match.arg(output)
  plot.data <- boostmtree.vimp.plot.data(
    x = x,
    show.interaction = show.interaction,
    show.time.effect = show.time.effect,
    use.percent = TRUE
  )
  if (identical(output, "data")) {
    return(invisible(plot.data))
  }
  if (identical(output, "pdf")) {
    if (is.null(file)) {
      file <- "boostmtree_vimp_plot.pdf"
    }
    pdf(file = file, width = 10, height = 7)
    on.exit(dev.off(), add = TRUE)
  }
  for (q in seq_along(plot.data)) {
    q.data <- plot.data[[q]]
    q.main <- pmax(q.data$main, 0)
    q.interaction <- if (!is.null(q.data$interaction)) pmax(q.data$interaction, 0) else NULL
    if (is.null(q.interaction)) {
      y.max <- max(q.main, na.rm = TRUE)
      if (!is.finite(y.max) || y.max <= 0) {
        y.max <- 1
      }
      barplot(
        height = q.main,
        names.arg = q.data$x.var.names,
        las = 2,
        col = col,
        ylim = c(0, y.max + eps),
        ylab = main,
        main = if (ncol(x$main) > 1L) paste(main, "-", q.data$q.name) else main,
        cex.names = cex.names,
        ...
      )
      if (!is.null(q.data$time.effect) && is.finite(q.data$time.effect)) {
        mtext(
          text = sprintf("Time effect: %.2f", q.data$time.effect),
          side = 3,
          line = 0.25,
          adj = 1,
          cex = 0.8
        )
      }
    } else {
      y.max <- max(c(q.main, q.interaction), na.rm = TRUE)
      if (!is.finite(y.max) || y.max <= 0) {
        y.max <- 1
      }
      barplot(
        height = q.main,
        names.arg = q.data$x.var.names,
        las = 2,
        col = col,
        ylim = c(-(y.max + eps), y.max + eps),
        ylab = main,
        yaxt = "n",
        main = if (ncol(x$main) > 1L) paste(main, "-", q.data$q.name) else main,
        cex.names = cex.names,
        ...
      )
      barplot(
        height = -q.interaction,
        col = col,
        add = TRUE,
        axes = FALSE,
        names.arg = FALSE,
        ...
      )
      abline(h = 0, lty = 1)
      if (!is.null(q.data$time.effect) && is.finite(q.data$time.effect)) {
        mtext(
          text = sprintf("Time effect: %.2f", q.data$time.effect),
          side = 3,
          line = 0.25,
          adj = 1,
          cex = 0.8
        )
      }
      mtext("Main effects", side = 3, line = -1.1, adj = 0, cex = 0.8)
      mtext("Time interactions", side = 1, line = 5, adj = 0, cex = 0.8)
      axis(2, at = pretty(c(-y.max, y.max)), labels = abs(pretty(c(-y.max, y.max))))
    }
  }
  invisible(x)
}
