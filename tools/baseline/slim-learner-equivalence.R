# Same-machine, bit-for-bit equivalence check for changes that only touch what
# a fitted boostmtree object stores (2.0.9: slimmer stored rfsrc trees and
# one-matrix-per-subject cross-validated coefficients).
#
# A storage change must not move any result. This script fits four models
# (cv.flag TRUE/FALSE x na.action "na.impute"/"na.omit", with a factor
# covariate and missing covariate values) and records every user-facing output:
# the fit's mu, err.rate, rmse, m.opt and gamma; predict() on the training data
# and on new data, with and without use.cv.flag; partial.plot() and
# marginal.plot() data; and vimp() of the grow and the predict object. It then
# compares two such fingerprints with identical().
#
# Like compare.R, it is meant to be run on ONE machine against two installed
# builds; across machines use reproducibility-check.R instead.
#
# usage:
#   Rscript tools/baseline/slim-learner-equivalence.R fingerprint <lib> <out.rds>
#   Rscript tools/baseline/slim-learner-equivalence.R compare <before.rds> <after.rds>
# <lib> is a library holding the boostmtree build to fingerprint.
# compare exits 0 when every component is identical, 1 otherwise.

args <- commandArgs(trailingOnly = TRUE)
usage <- "usage: slim-learner-equivalence.R fingerprint <lib> <out.rds> | compare <before.rds> <after.rds>"
if (length(args) != 3L || !args[1] %in% c("fingerprint", "compare")) {
  stop(usage)
}

if (args[1] == "compare") {
  before <- readRDS(args[2])
  after <- readRDS(args[3])
  if (length(before) == 0L || !identical(names(before), names(after))) {
    cat("FAIL: the fingerprints hold different components\n")
    quit(status = 1L)
  }
  failed <- character(0)
  for (nm in names(before)) {
    ok <- identical(before[[nm]], after[[nm]])
    if (is.null(before[[nm]])) {
      ok <- FALSE
    }
    cat(sprintf("  %-40s %s\n", nm, if (ok) "identical" else "DIFFERS"))
    if (!ok) failed <- c(failed, nm)
  }
  if (length(failed) > 0L) {
    cat("FAIL:", length(failed), "of", length(before), "components differ\n")
    quit(status = 1L)
  }
  cat("PASS: all", length(before), "components identical\n")
  quit(status = 0L)
}

suppressMessages(library(boostmtree, lib.loc = args[2]))
cat("boostmtree", as.character(packageVersion("boostmtree", lib.loc = args[2])), "from", args[2], "\n")

# Drop what legitimately differs between two runs of the same code: matched
# calls, and the stored trees themselves (the change under test).
scrub <- function(x) {
  if (is.list(x)) {
    keep <- setdiff(names(x), c("call", "base.learner", "base.learner.xvar"))
    if (!is.null(names(x))) {
      cls <- class(x)
      x <- lapply(unclass(x)[keep], scrub)
      class(x) <- cls
    } else {
      x <- lapply(x, scrub)
    }
  }
  x
}

# The per-subject cross-validated coefficients, read through the package's
# own parser so the comparison is about values, not storage layout.
gamma.i.values <- function(gamma.i.list) {
  parse <- getFromNamespace("boostmtree.predict.normalize.gamma.step", "boostmtree")
  lapply(gamma.i.list, function(by.m) {
    lapply(by.m, function(by.i) {
      lapply(by.i, function(step) {
        if (is.null(step)) NULL else lapply(parse(step), function(v) as.numeric(v))
      })
    })
  })
}

set.seed(20261113)
d <- simLong(n = 60, n.time = 3, rho = 0.8, model = 2, family = "continuous", q = 4)$data.list
x <- d$features
x$grp <- factor(sample(c("a", "b", "c"), nrow(x), replace = TRUE))
subject.rows <- split(seq_len(nrow(x)), d$id)
# Subject-level covariates must be constant within subject: knock out whole
# subjects' values so NA handling is exercised the way it is in practice.
for (s in c(3, 11, 29)) x$x1[subject.rows[[s]]] <- NA
for (s in c(7, 41)) x$grp[subject.rows[[s]]] <- NA

new.ids <- unique(d$id)[1:15]
new.rows <- d$id %in% new.ids

out <- list()
for (cv in c(TRUE, FALSE)) {
  for (na in c("na.impute", "na.omit")) {
    tag <- paste0("cv=", cv, ",", na)
    set.seed(11)
    fit <- boostmtree(x, d$time, d$id, d$y, family = "continuous", M = 12,
                      cv.flag = cv, na.action = na, mtry = 3, verbose = FALSE)
    out[[paste(tag, "fit")]] <- scrub(fit[c("mu", "err.rate", "rmse", "m.opt", "gamma",
                                           "lambda", "phi", "rho", "membership")])
    if (cv) {
      out[[paste(tag, "gamma.i")]] <- gamma.i.values(fit$gamma.i.list)
    }
    out[[paste(tag, "predict.train")]] <- scrub(predict(fit))
    if (cv) {
      # This predict object carries the per-subject coefficients as `gamma`;
      # compare them by value, like gamma.i above.
      pred.cv <- scrub(predict(fit, use.cv.flag = TRUE))
      pred.cv$gamma <- gamma.i.values(pred.cv$gamma)
      out[[paste(tag, "predict.train.cv")]] <- pred.cv
    }
    pred.new <- predict(fit, x = x[new.rows, ], tm = d$time[new.rows],
                        id = d$id[new.rows], y = d$y[new.rows])
    out[[paste(tag, "predict.new")]] <- scrub(pred.new)
    set.seed(5)
    out[[paste(tag, "partial.plot")]] <- scrub(partial.plot(
      fit, x.var.names = c("x2", "grp"), time.points = c(0.25, 0.75), n.points = 4,
      output = "data", verbose = FALSE
    ))
    set.seed(5)
    out[[paste(tag, "marginal.plot")]] <- scrub(marginal.plot(
      fit, x.var.names = c("x2", "grp"), time.points = c(0.25, 0.75),
      output = "data", verbose = FALSE
    ))
    set.seed(5)
    out[[paste(tag, "vimp.predict")]] <- vimp.boostmtree(pred.new)
    if (cv) {
      set.seed(5)
      out[[paste(tag, "vimp.grow")]] <- vimp.boostmtree(fit)
      set.seed(5)
      out[[paste(tag, "vimp.grow.joint")]] <- vimp.boostmtree(fit, joint = TRUE)
    }
  }
}
saveRDS(out, args[3])
cat("wrote", length(out), "components to", args[3], "\n")
