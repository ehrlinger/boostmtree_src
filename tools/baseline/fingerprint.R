# Deterministic fingerprint of the cv.flag = FALSE fitting path.
# Used to prove the cv.flag fix introduces no collateral change.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript fingerprint.R <output.rds>")

suppressMessages(library(boostmtree))

set.seed(20260902)
d <- simLong(
  n = 60, n.time = 5, rho = 0.8, model = 2,
  family = "continuous", q = 0
)$data.list

fit <- boostmtree(
  d$features, d$time, d$id, d$y,
  family = "continuous", M = 50, cv.flag = FALSE
)

saveRDS(
  list(mu = fit$mu, lambda = fit$lambda, gamma = fit$gamma,
       phi = fit$phi, rho = fit$rho),
  args[1]
)
cat("fingerprint written to", args[1], "\n")
