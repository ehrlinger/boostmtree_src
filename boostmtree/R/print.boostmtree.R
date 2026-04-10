print.boostmtree <- function(x, ...) {
  display <- boostmtree.resolve.display.object(x)
  model <- display$model
  univariate <- isTRUE(model$univariate)
  cat("model                       :", class(model)[3], "\n")
  cat("fitting mode                :", display$kind, "\n")
  cat("family                      :", model$family, "\n")
  if (model$family %in% c("binary", "nominal", "ordinal")) {
    cat("number of levels            :", model$q.total, "\n")
  }
  if (!is.null(model$ntree) && model$ntree > 1L) {
    cat("ntree                       :", model$ntree, "\n")
  }
  cat("number of terminal nodes    :", model$k, "\n")
  cat("regularization parameter    :", model$nu[1], "\n")
  cat("sample size                 :", nrow(display$x), "\n")
  cat("number of variables         :", ncol(display$x), "\n")
  if (!univariate) {
    cat("number of unique time points:", length(model$time.unique), "\n")
    cat("avg. number of time points  :", round(mean(vapply(display$time, length, integer(1))), 2), "\n")
    cat("B-spline dimension          :", ncol(model$x.tm), "\n")
    cat("penalization order          :", model$pen.ord, "\n")
  } else {
    cat("univariate family           :", TRUE, "\n")
  }
  cat("boosting iterations         :", model$M, "\n")
  if (!is.null(display$err.rate) && !is.null(display$m.opt)) {
    optimized.rho <- if (!is.null(model$rho) && !is.null(display$m.opt)) {
      vapply(seq_len(model$n.q), function(q) {
        if (is.matrix(model$rho)) model$rho[display$m.opt[q], q] else model$rho[display$m.opt[q]]
      }, numeric(1))
    } else {
      NULL
    }
    optimized.phi <- if (!is.null(model$phi) && !is.null(display$m.opt)) {
      vapply(seq_len(model$n.q), function(q) {
        if (is.matrix(model$phi)) model$phi[display$m.opt[q], q] else model$phi[display$m.opt[q]]
      }, numeric(1))
    } else {
      NULL
    }
    cat("optimized number iterations :", display$m.opt, "\n")
    if (!univariate && !is.null(optimized.rho)) {
      cat("optimized rho               :", round(optimized.rho, 4), "\n")
    }
    if (!univariate && !is.null(optimized.phi)) {
      cat("optimized phi               :", round(optimized.phi, 4), "\n")
    }
    metric.label <- if (identical(display$kind, "grow")) {
      "OOB cv RMSE"
    } else {
      "test set RMSE"
    }
    cat(sprintf("%-28s:", metric.label), round(display$rmse, 4), "\n")
  }
  invisible(x)
}
