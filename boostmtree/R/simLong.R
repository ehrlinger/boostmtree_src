simLong.random.correlation <- function(n, rho) {
  raw <- matrix(stats::runif(n * n, min = -abs(rho), max = abs(rho)), nrow = n)
  raw <- (raw + t(raw)) / 2
  diag(raw) <- 1
  eig <- eigen(raw, symmetric = TRUE)
  eig$values[eig$values < 1e-6] <- 1e-6
  sigma <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  scale.inv <- diag(1 / sqrt(diag(sigma)), nrow = n)
  scale.inv %*% sigma %*% scale.inv
}
simLong <- function(
  n = 100,
  n.test = 0,
  n.time = 5,
  rho = 0.8,
  cor.type = c("cor.comp.sym", "cor.ar1", "cor.symm", "iid"),
  model = c(0, 1, 2, 3),
  family = c("continuous", "binary"),
  phi = 1,
  q = 0,
  ...
) {
  family <- match.arg(family)
  cor.type <- match.arg(cor.type)
  model <- as.integer(match.arg(as.character(model), as.character(0:3)))
  if (length(n.time) != 1L || is.na(n.time) || n.time < 1) {
    stop("`n.time` must be a single positive integer.")
  }
  n.time <- as.integer(n.time)
  if (length(n.test) != 1L || is.na(n.test) || n.test < 0) {
    stop("`n.test` must be a single non-negative integer.")
  }
  n.test <- as.integer(n.test)
  if (length(phi) != 1L || is.na(phi) || phi <= 0) {
    stop("`phi` must be a single positive numeric value.")
  }
  phi <- as.numeric(phi)
  if (length(rho) != 1L || is.na(rho) || abs(rho) >= 1) {
    stop("`rho` must be a single numeric value strictly inside (-1, 1).")
  }
  rho <- as.numeric(rho)
  if (length(q) != 1L || is.na(q) || q < 0) {
    stop("`q` must be a single non-negative integer.")
  }
  q <- as.integer(q)
  data <- data.frame(do.call("rbind", lapply(seq_len(n + n.test), function(i) {
    n.i <- round(stats::runif(1, min = 1, max = 3 * n.time))
    corr.matrix <- switch(
      cor.type,
      "cor.comp.sym" = {
        out <- matrix(rho, nrow = n.i, ncol = n.i)
        diag(out) <- 1
        out
      },
      "cor.ar1" = {
        out <- diag(rep(1, n.i))
        if (n.i > 1L) {
          for (j in seq_len(n.i - 1L)) {
            out[j, (j + 1L):n.i] <- rho^(seq_len(n.i - j))
          }
          out[lower.tri(out)] <- t(out)[lower.tri(out)]
        }
        out
      },
      "cor.symm" = simLong.random.correlation(n.i, rho),
      "iid" = diag(rep(1, n.i))
    )
    eps <- sqrt(phi) * t(chol(corr.matrix)) %*% stats::rnorm(n.i)
    x1 <- stats::rnorm(1)
    x2 <- stats::runif(1, 1, 2)
    x3 <- stats::runif(1, 1, 3)
    x4 <- stats::rnorm(1)
    x <- c(x1, x2, x3, x4)
    if (q > 0L) {
      x <- c(x, stats::rnorm(q))
    }
    tm <- sample((seq_len(3 * n.time)) / n.time, size = n.i, replace = TRUE)
    linear.predictor <- switch(
      as.character(model),
      "0" = 1.5 + 2.5 * x1 - 1.2 * x3 - 0.6 * x4 + eps,
      "1" = 1.5 + 2.5 * x1 - 1.2 * x3 - 0.2 * x4 - 0.65 * tm * x2 + eps,
      "2" = 1.5 + 2.5 * x1 - 1.2 * x3 - 0.2 * x4 - 0.65 * (tm^2) * (x2^2) + eps,
      "3" = 1.5 + 2.5 * x1 - 1.2 * x3 - 0.2 * exp(x4) - 0.65 * (tm^2) * (x2^2) * x3 + eps
    )
    mu <- boostmtree.get.mu(linear.predictor = linear.predictor, family = family)
    y <- if (family == "continuous") {
      mu
    } else {
      stats::rbinom(n = n.i, size = 1, prob = mu)
    }
    cbind(
      matrix(x, nrow = n.i, ncol = length(x), byrow = TRUE),
      tm,
      rep(i, n.i),
      y
    )
  })))
  n.feature <- q + 4L
  colnames(data) <- c(paste0("x", seq_len(n.feature)), "time", "id", "y")
  data.list <- list(
    features = data[, seq_len(n.feature), drop = FALSE],
    time = data$time,
    id = data$id,
    y = data$y
  )
  formula.true <- switch(
    as.character(model),
    "0" = "y ~ g(x1 + x3 + x4)",
    "1" = "y ~ g(x1 + x3 + x4 + I(time * x2))",
    "2" = "y ~ g(x1 + x3 + x4 + I(time^2 * x2^2))",
    "3" = "y ~ g(x1 + x3 + exp(x4) + I(time^2 * x2^2 * x3))"
  )
  train.index <- which(data$id <= n)
  invisible(list(
    data.list = data.list,
    data = data,
    train.index = train.index,
    formula.true = formula.true
  ))
}
