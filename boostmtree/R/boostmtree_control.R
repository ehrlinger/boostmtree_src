boostmtree.control <- function(
  ntree = 1,
  bootstrap = "by.root",
  bootstrap.fraction = 0.632,
  sample.matrix = NULL,
  nsplit = NULL,
  samptype = "swor",
  xvar.wt = NULL,
  case.wt = NULL,
  seed = NULL,
  cv.lambda = FALSE,
  cv.rho = TRUE
) {
  structure(
    list(
      ntree = ntree,
      bootstrap = bootstrap,
      bootstrap.fraction = bootstrap.fraction,
      sample.matrix = sample.matrix,
      nsplit = nsplit,
      samptype = samptype,
      xvar.wt = xvar.wt,
      case.wt = case.wt,
      seed = seed,
      cv.lambda = cv.lambda,
      cv.rho = cv.rho
    ),
    class = c("boostmtree.control", "list")
  )
}
boostmtree.normalize.control <- function(control) {
  default.control <- boostmtree.control()
  if (is.null(control)) {
    control <- default.control
  }
  if (!inherits(control, "boostmtree.control")) {
    if (!is.list(control)) {
      stop("`control` must be NULL, a list, or a `boostmtree.control` object.")
    }
    control <- modifyList(default.control, control)
    class(control) <- class(default.control)
  } else {
    control <- modifyList(default.control, control)
    class(control) <- class(default.control)
  }
  if (length(control$ntree) != 1L || is.na(control$ntree)) {
    stop("`control$ntree` must be a single positive integer.")
  }
  control$ntree <- as.integer(control$ntree)
  if (control$ntree < 1L) {
    stop("`control$ntree` must be at least 1.")
  }
  if (control$ntree != 1L) {
    stop(
      "This first-pass refactor supports `control$ntree = 1` only. ",
      "The forest path is intentionally deferred."
    )
  }
  if (length(control$bootstrap) != 1L || !nzchar(control$bootstrap)) {
    stop("`control$bootstrap` must be a non-empty character scalar.")
  }
  control$bootstrap <- match.arg(
    as.character(control$bootstrap),
    c("by.root", "by.user", "none")
  )
  if (length(control$bootstrap.fraction) != 1L || is.na(control$bootstrap.fraction)) {
    stop("`control$bootstrap.fraction` must be a single numeric value.")
  }
  control$bootstrap.fraction <- as.numeric(control$bootstrap.fraction)
  if (control$bootstrap.fraction <= 0 || control$bootstrap.fraction > 1) {
    stop("`control$bootstrap.fraction` must lie in (0, 1].")
  }
  if (!is.null(control$sample.matrix) && !is.matrix(control$sample.matrix)) {
    stop("`control$sample.matrix` must be NULL or a matrix.")
  }
  if (!is.null(control$nsplit)) {
    if (length(control$nsplit) != 1L || is.na(control$nsplit)) {
      stop("`control$nsplit` must be NULL or a single positive integer.")
    }
    control$nsplit <- as.integer(control$nsplit)
    if (control$nsplit < 1L) {
      stop("`control$nsplit` must be at least 1.")
    }
  }
  if (length(control$samptype) != 1L || !nzchar(control$samptype)) {
    stop("`control$samptype` must be a non-empty character scalar.")
  }
  control$samptype <- as.character(control$samptype)
  if (!is.null(control$seed)) {
    if (length(control$seed) != 1L || is.na(control$seed)) {
      stop("`control$seed` must be NULL or a single integer.")
    }
    control$seed <- as.integer(control$seed)
  }
  control$cv.lambda <- isTRUE(control$cv.lambda)
  control$cv.rho <- isTRUE(control$cv.rho)
  if (identical(control$bootstrap, "none") && !is.null(control$sample.matrix)) {
    warning(
      "`control$sample.matrix` is ignored when `control$bootstrap = \"none\"`."
    )
    control$sample.matrix <- NULL
  }
  control
}
boostmtree.oob.count.from.sample.matrix <- function(sample.matrix) {
  if (is.null(sample.matrix)) {
    return(integer(0))
  }
  colSums(sample.matrix == 0L)
}
boostmtree.enforce.oob.control <- function(control, cv.flag, n, M) {
  if (!isTRUE(cv.flag)) {
    return(control)
  }
  if (identical(control$bootstrap, "none")) {
    warning(
      "`cv.flag = TRUE` uses out-of-bag error. Replacing ",
      "`control$bootstrap = \"none\"` with `control$bootstrap = \"by.root\"` ",
      "for this fit."
    )
    control$bootstrap <- "by.root"
    control$sample.matrix <- NULL
    return(control)
  }
  if (identical(control$bootstrap, "by.user")) {
    if (is.null(control$sample.matrix)) {
      if (control$bootstrap.fraction >= 1) {
        warning(
          "`cv.flag = TRUE` requires out-of-bag subjects. Replacing ",
          "`control$bootstrap = \"by.user\"` and ",
          "`control$bootstrap.fraction = 1` with `control$bootstrap = \"by.root\"` ",
          "for this fit."
        )
        control$bootstrap <- "by.root"
      }
      return(control)
    }
    sample.matrix <- boostmtree.resolve.sample.matrix(control = control, n = n, M = M)
    oob.count <- boostmtree.oob.count.from.sample.matrix(sample.matrix)
    if (any(oob.count == 0L)) {
      stop(
        "`cv.flag = TRUE` requires at least one out-of-bag subject at each boosting iteration. ",
        "`control$sample.matrix` has no out-of-bag subjects for iteration(s): ",
        paste(which(oob.count == 0L), collapse = ", "),
        ". Use the default `control = boostmtree.control()` or provide a sample matrix ",
        "with at least one zero in every column."
      )
    }
  }
  control
}
boostmtree.build.sample.matrix <- function(n, M, bootstrap.fraction) {
  sample.matrix <- matrix(NA_integer_, nrow = n, ncol = M)
  for (m in seq_len(M)) {
    sample.index <- sample.int(n, size = floor(bootstrap.fraction * n), replace = FALSE)
    sample.index <- sort(c(
      sample.index,
      sample(sample.index, size = n - length(sample.index), replace = TRUE)
    ))
    sample.matrix[, m] <- vapply(
      seq_len(n),
      function(i) sum(sample.index == i),
      integer(1)
    )
  }
  sample.matrix
}
boostmtree.resolve.sample.matrix <- function(control, n, M) {
  if (!identical(control$bootstrap, "by.user")) {
    return(NULL)
  }
  if (is.null(control$sample.matrix)) {
    return(boostmtree.build.sample.matrix(
      n = n,
      M = M,
      bootstrap.fraction = control$bootstrap.fraction
    ))
  }
  sample.matrix <- control$sample.matrix
  if (!identical(dim(sample.matrix), c(n, M))) {
    stop(
      "`control$sample.matrix` must have dimensions ",
      n, " x ", M, "."
    )
  }
  sample.matrix
}
