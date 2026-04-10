boostmtree.split.by.subject <- function(values, id, id.unique) {
  unname(split(values, factor(id, levels = id.unique)))
}
boostmtree.remove.missing.subjects <- function(x.subject, na.action) {
  if (ncol(x.subject) == 0L) {
    stop("`x` must contain at least one subject-level covariate.")
  }
  keep.columns <- which(vapply(
    seq_len(ncol(x.subject)),
    function(j) !all(is.na(x.subject[[j]])),
    logical(1)
  ))
  if (length(keep.columns) == 0L) {
    stop("No usable subject-level covariates remain after removing all-missing columns.")
  }
  x.reduced <- x.subject[, keep.columns, drop = FALSE]
  row.all.na <- vapply(
    seq_len(nrow(x.reduced)),
    function(i) all(is.na(x.reduced[i, , drop = TRUE])),
    logical(1)
  )
  row.any.na <- vapply(
    seq_len(nrow(x.reduced)),
    function(i) any(is.na(x.reduced[i, , drop = TRUE])),
    logical(1)
  )
  remove.row <- if (identical(na.action, "na.omit")) row.any.na else row.all.na
  keep.subject <- which(!remove.row)
  if (length(keep.subject) == 0L) {
    stop("No subjects remain after applying the missing-data rule.")
  }
  list(
    x.subject = x.reduced[keep.subject, , drop = FALSE],
    keep.columns = keep.columns,
    keep.subject = keep.subject
  )
}
boostmtree.preprocess.data <- function(x, tm = NULL, id = NULL, y, na.action = c("na.omit", "na.impute")) {
  na.action <- match.arg(na.action)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  n.obs <- nrow(x)
  if (length(y) != n.obs) {
    stop("`y` must have length equal to `nrow(x)`.")
  }
  if (is.null(tm)) {
    id <- seq_len(n.obs)
    tm <- rep(0, n.obs)
  } else {
    if (is.null(id)) {
      stop("`id` must be supplied when `tm` is supplied.")
    }
    if (length(tm) != n.obs) {
      stop("`tm` must have length equal to `nrow(x)`.")
    }
    if (length(id) != n.obs) {
      stop("`id` must have length equal to `nrow(x)`.")
    }
  }
  if (anyNA(id) || anyNA(y) || anyNA(tm)) {
    stop("Missing values in `id`, `y`, or `tm` are not allowed.")
  }
  id.unique <- sort(unique(id))
  univariate <- length(id.unique) == length(id)
  if (univariate) {
    tm <- rep(0, n.obs)
  }
  row.order <- order(match(id, id.unique), seq_along(id))
  x.long <- x[row.order, , drop = FALSE]
  tm <- tm[row.order]
  id <- id[row.order]
  y <- y[row.order]
  id.unique <- sort(unique(id))
  subject.start <- match(id.unique, id)
  x.subject.raw <- x.long[subject.start, , drop = FALSE]
  missing.info <- boostmtree.remove.missing.subjects(x.subject.raw, na.action = na.action)
  kept.ids <- id.unique[missing.info$keep.subject]
  keep.long <- id %in% kept.ids
  x.long <- x.long[keep.long, missing.info$keep.columns, drop = FALSE]
  tm <- tm[keep.long]
  id <- id[keep.long]
  y <- y[keep.long]
  id.unique <- kept.ids
  subject.start <- match(id.unique, id)
  x.subject <- x.long[subject.start, , drop = FALSE]
  n.subject <- length(id.unique)
  ni <- unname(tabulate(match(id, id.unique), nbins = n.subject))
  time.by.subject <- boostmtree.split.by.subject(tm, id, id.unique)
  y.by.subject <- boostmtree.split.by.subject(y, id, id.unique)
  list(
    x.subject = x.subject,
    x.long = x.long,
    time = tm,
    id = id,
    y = y,
    id.unique = id.unique,
    n = n.subject,
    ni = ni,
    time.by.subject = time.by.subject,
    y.by.subject = y.by.subject,
    x.var.names = colnames(x.subject),
    univariate = univariate
  )
}
boostmtree.build.time.design <- function(time.by.subject, d, n.knots) {
  time.unique <- sort(unique(unlist(time.by.subject)))
  if (n.knots < 0) {
    warning("B-splines require a non-negative number of knots; using an intercept-only time design.")
    d <- 0
  }
  if (d >= 1) {
    if (length(time.unique) > 1L) {
      bs.time <- bs(time.unique, df = n.knots + d, degree = d)
      x.tm <- cbind(1, bs.time)
      attr(x.tm, "knots") <- attr(bs.time, "knots")
      attr(x.tm, "Boundary.knots") <- attr(bs.time, "Boundary.knots")
    } else {
      x.tm <- cbind(1, cbind(time.unique))
    }
  } else {
    x.tm <- cbind(rep(1, length(time.unique)))
  }
  time.index.by.subject <- lapply(time.by.subject, function(time.i) {
    match(time.i, time.unique)
  })
  time.design <- lapply(time.index.by.subject, function(index.i) {
    x.tm[index.i, , drop = FALSE]
  })
  list(
    time.unique = time.unique,
    x.tm = x.tm,
    time.design = time.design,
    d = d,
    df.time.design = ncol(x.tm)
  )
}
