# Cross-machine reproducibility check for the boostmtree CCF fork.
#
# WHY THIS EXISTS
# ---------------
# boostmtree fits by asking randomForestSRC for terminal-node memberships
# (predict.rfsrc(..., ptn.count = k)$ptn.membership). DESCRIPTION only requires
# randomForestSRC (>= 3.5.0), but every verification run for this fork was done
# against randomForestSRC 3.6.2 on macOS arm64 under R 4.6.1. A different
# randomForestSRC could partition the covariate space differently, which changes
# gamma, which changes the fitted values -- silently, with no error raised.
# This script re-runs the reference fit locally and tells you whether your
# environment reproduces the reference numbers before you trust them.
#
# It is the tolerant sibling of compare.R. compare.R uses identical() and is the
# right tool for proving a code change moved nothing ON ONE MACHINE. Run ACROSS
# machines, identical() is the wrong test: a different BLAS/LAPACK, a different
# CPU, or FMA contraction will legitimately move the last bits of a correct
# result and identical() would cry wolf. compare.R is left exactly as it was.
#
# usage: Rscript tools/baseline/reproducibility-check.R [baseline.rds]
#   default baseline: tools/baseline/baseline-cvflag-false.rds (next to this file)
# exit 0 = your environment reproduces the reference fit within tolerance
# exit 1 = it does not; see the failure message.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("usage: Rscript reproducibility-check.R [baseline.rds]")
}

## ---------------------------------------------------------------------------
## Tolerance
## ---------------------------------------------------------------------------
# 1e-8 relative, applied through all.equal() (mean relative difference).
#
# Rationale. IEEE doubles carry ~2.2e-16 of relative resolution. A 50-stage
# boosted fit with a penalised-spline solve accumulates rounding over many
# thousands of flops, so two *correct* runs on different BLAS/LAPACK builds or
# different CPUs (fused multiply-add, different reduction orders, different
# vector widths) routinely disagree in the last 3-6 significant digits. 1e-8 is
# roughly sqrt(.Machine$double.eps) -- half of machine precision -- and is the
# same order as all.equal()'s own default (1.5e-8). It is therefore loose enough
# never to fire on honest floating-point noise.
#
# It is also far tighter than any failure this check is meant to catch. A
# different randomForestSRC version producing different terminal-node
# memberships does not perturb gamma in the 8th digit; it relocates whole nodes
# and moves fitted values by O(1) relative. (For scale: the cv.flag bug this
# fork fixes inflated terminal-node coefficients by ~5x.) So a real divergence
# overshoots 1e-8 by seven or more orders of magnitude, and there is no useful
# grey zone between "last-bit noise" and "different model".
TOL <- 1e-8

## ---------------------------------------------------------------------------
## Environment report -- print this FIRST so it is captured even on a crash.
## ---------------------------------------------------------------------------
si <- sessionInfo()
cat("boostmtree cross-machine reproducibility check\n")
cat("=============================================\n\n")
cat("Local environment\n")
cat("  R version       :", R.version.string, "\n")
cat("  Platform        :", R.version$platform, "\n")
cat("  Running under   :", si$running, "\n")
pkg.version <- function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) "<not installed>" else v
}
cat("  boostmtree      :", pkg.version("boostmtree"), "\n")
cat("  randomForestSRC :", pkg.version("randomForestSRC"), "\n")
cat("  nlme            :", pkg.version("nlme"), "\n")
cat("  BLAS            :", if (is.null(si$BLAS)) "<unknown>" else si$BLAS, "\n")
cat("  LAPACK          :", if (is.null(si$LAPACK)) "<unknown>" else si$LAPACK, "\n\n")
cat("Reference environment (where this fork was verified)\n")
cat("  R 4.6.1 / aarch64-apple-darwin23 / randomForestSRC 3.6.2 / nlme 3.1-169\n\n")

## ---------------------------------------------------------------------------
## Locate the baseline
## ---------------------------------------------------------------------------
script.dir <- {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
}
baseline.path <- if (length(args) == 1L) {
  args[1]
} else {
  file.path(script.dir, "baseline-cvflag-false.rds")
}
if (!file.exists(baseline.path)) {
  cat("ERROR: baseline fingerprint not found:", baseline.path, "\n")
  quit(status = 1L)
}
cat("Baseline:", baseline.path, "\n\n")
baseline <- readRDS(baseline.path)

## ---------------------------------------------------------------------------
## Guard against a vacuous pass (same contract as compare.R)
## ---------------------------------------------------------------------------
components <- names(baseline)
if (length(components) == 0L) {
  cat("ERROR: baseline fingerprint has no components; nothing was compared.\n")
  quit(status = 1L)
}
null.components <- components[vapply(components,
                                     function(nm) is.null(baseline[[nm]]),
                                     logical(1))]
if (length(null.components) > 0L) {
  cat("ERROR: baseline component(s) NULL, comparison would be vacuous:",
      paste(null.components, collapse = ", "), "\n")
  quit(status = 1L)
}

## ---------------------------------------------------------------------------
## Reproduce the reference fit.
## Simulation and fit settings are copied verbatim from fingerprint.R so the
## two fingerprints are directly comparable. Do not change them here without
## regenerating the committed baseline.
## ---------------------------------------------------------------------------
suppressMessages(library(boostmtree))

set.seed(20260902)
d <- simLong(
  n = 60, n.time = 5, rho = 0.8, model = 2,
  family = "continuous", q = 0
)$data.list

cat("Refitting the reference model (n = 60, M = 50, cv.flag = FALSE) ...\n")
invisible(utils::capture.output(
  fit <- boostmtree(
    d$features, d$time, d$id, d$y,
    family = "continuous", M = 50, cv.flag = FALSE
  )
))

current <- list(mu = fit$mu, lambda = fit$lambda, gamma = fit$gamma,
                phi = fit$phi, rho = fit$rho)
cat("done.\n\n")

## ---------------------------------------------------------------------------
## Compare, component by component
## ---------------------------------------------------------------------------
# The fingerprint components are nested lists of numeric vectors / matrices.
# Flattening to a plain numeric vector lets us take real max-abs and max-rel
# differences; the length check below catches any structural divergence that
# flattening would otherwise hide.
flatten.numeric <- function(x) {
  suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
}

failed <- character(0)
rows <- list()

for (nm in components) {
  b <- flatten.numeric(baseline[[nm]])
  p <- flatten.numeric(current[[nm]])

  if (length(b) != length(p)) {
    failed <- c(failed, nm)
    rows[[nm]] <- sprintf("  %-8s STRUCTURE MISMATCH  baseline has %d values, local fit has %d",
                          nm, length(b), length(p))
    next
  }
  if (length(b) == 0L) {
    failed <- c(failed, nm)
    rows[[nm]] <- sprintf("  %-8s EMPTY               component holds no numeric values", nm)
    next
  }
  if (anyNA(b) != anyNA(p) || (anyNA(p) && !identical(is.na(b), is.na(p)))) {
    failed <- c(failed, nm)
    rows[[nm]] <- sprintf("  %-8s NA PATTERN MISMATCH baseline NAs = %d, local NAs = %d",
                          nm, sum(is.na(b)), sum(is.na(p)))
    next
  }

  finite <- is.finite(b) & is.finite(p)
  abs.diff <- abs(p[finite] - b[finite])
  max.abs <- if (length(abs.diff)) max(abs.diff) else 0
  nz <- finite & b != 0
  max.rel <- if (any(nz)) max(abs(p[nz] - b[nz]) / abs(b[nz])) else NA_real_

  cmp <- all.equal(b, p, tolerance = TOL)
  ok <- isTRUE(cmp)
  if (!ok) failed <- c(failed, nm)

  rows[[nm]] <- sprintf("  %-8s %-5s n = %-6d max|diff| = %-11.4g max rel diff = %-11.4g%s",
                        nm, if (ok) "OK" else "FAIL", length(b), max.abs,
                        if (is.na(max.rel)) NA_real_ else max.rel,
                        if (ok) "" else paste0("\n           all.equal(): ",
                                               paste(cmp, collapse = "; ")))
}

cat("Per-component comparison (tolerance =", format(TOL), "relative)\n")
for (nm in components) cat(rows[[nm]], "\n")
cat("\n")

## ---------------------------------------------------------------------------
## Verdict
## ---------------------------------------------------------------------------
if (length(failed) == 0L) {
  cat("PASS: this environment reproduces the reference cv.flag = FALSE fit\n")
  cat("      within", format(TOL), "relative tolerance across", length(components),
      "components:", paste(components, collapse = ", "), "\n")
  quit(status = 0L)
}

cat("FAIL: this environment does NOT reproduce the reference fit.\n")
cat("      Divergent component(s):", paste(failed, collapse = ", "), "\n\n")
cat("What this means\n")
cat("  Your machine produces different fitted values from the environment this\n")
cat("  fork was verified in. Differences of this size are not floating-point\n")
cat("  noise, so mu, gamma, predict(), partial.plot() and vimp() output from\n")
cat("  this installation may not match the reference results.\n\n")
cat("What to do\n")
cat("  1. Compare the environment block above with the reference environment.\n")
cat("     The usual cause is a different randomForestSRC: boostmtree reads\n")
cat("     terminal-node memberships from predict.rfsrc(..., ptn.count = k), and\n")
cat("     a different version can partition the data differently.\n")
cat("  2. Pin randomForestSRC to the reference version before trusting numbers:\n")
cat("       remotes::install_version(\"randomForestSRC\", version = \"3.6.2\")\n")
cat("  3. Re-run this script. If it still fails with randomForestSRC pinned,\n")
cat("     do not use this build for the analysis -- report the output above to\n")
cat("     the package maintainer.\n")
quit(status = 1L)
