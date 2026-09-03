# Post-install smoke test for the boostmtree CCF fork.
#
# Run against an INSTALLED boostmtree (not the source tree) to confirm the
# installed build actually exercises the fixed code path. Deliberately small so
# it finishes in seconds on a CI runner.
#
# The assertion mirrors the style of
# boostmtree/tests/testthat/test-cv-flag-gamma-decay.R: it does not pin stored
# numerics (those are platform-dependent), it checks that the in-sample boost
# converged, so it stays valid on any OS. If the cv.flag = TRUE in-sample mu
# freeze were present, the linear predictor would accumulate M unshrunk
# increments and diverge well outside the observed response range.
#
# usage: Rscript tools/smoke/install-smoke.R
# exit 0 = installed build behaves; non-zero = it does not.

suppressMessages(library(boostmtree))

cat("boostmtree", as.character(utils::packageVersion("boostmtree")),
    "installed at", dirname(system.file(package = "boostmtree")), "\n")
cat("randomForestSRC", as.character(utils::packageVersion("randomForestSRC")), "\n")
cat(R.version.string, R.version$platform, "\n\n")

set.seed(202609)
d <- simLong(
  n = 50, n.time = 5, rho = 0.8, model = 2,
  family = "continuous", q = 0
)$data.list

fit <- boostmtree(
  d$features, d$time, d$id, d$y,
  family = "continuous", M = 30, cv.flag = TRUE, verbose = FALSE
)

# predict(use.cv.flag = FALSE) reads the in-sample path -- the one the cv.flag
# fix repaired, and the one partial.plot() and vimp() also read.
predicted.mu <- unlist(predict(fit, use.cv.flag = FALSE)$mu)

fail <- function(...) { cat("SMOKE TEST FAILED:", ..., "\n"); quit(status = 1L) }

if (length(predicted.mu) == 0L) fail("predict() returned no fitted values")
if (!all(is.finite(predicted.mu))) fail("predict() returned non-finite fitted values")

observed <- range(d$y)
width <- diff(observed)
lower <- observed[1] - 0.25 * width
upper <- observed[2] + 0.25 * width
predicted <- range(predicted.mu)

cat(sprintf("observed  y range : [%.4f, %.4f]\n", observed[1], observed[2]))
cat(sprintf("sane band         : [%.4f, %.4f]\n", lower, upper))
cat(sprintf("predicted mu range: [%.4f, %.4f]  (n = %d)\n",
            predicted[1], predicted[2], length(predicted.mu)))

if (predicted[1] < lower || predicted[2] > upper) {
  fail("predicted mu escapes the sane band around the observed response range;",
       "the cv.flag = TRUE in-sample path is diverging in this build")
}

cat("\nSMOKE TEST PASSED: installed build fits, predicts, and stays in range.\n")
