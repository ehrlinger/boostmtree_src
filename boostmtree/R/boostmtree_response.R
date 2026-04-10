boostmtree.match.family <- function(family) {
  match.arg(family, c("continuous", "binary", "nominal", "ordinal"))
}
boostmtree.prepare.response <- function(y, family, y.reference = NULL) {
  family <- boostmtree.match.family(family)
  if (family == "continuous") {
    if (!is.numeric(y)) {
      stop("`y` must be numeric for the continuous family.")
    }
    return(list(
      y.input = y,
      y.levels = NA,
      q.total = 1L,
      n.q = 1L,
      q.set.index = NA_integer_,
      q.set = NA,
      y.reference = NULL,
      y.by.q.vector = list(y)
    ))
  }
  if (is.factor(y)) {
    y.factor <- droplevels(y)
    y.levels <- levels(y.factor)
    y.value <- as.character(y.factor)
    if (!is.null(y.reference) && !is.character(y.reference)) {
      y.reference <- as.character(y.reference)
    }
  } else {
    y.value <- y
    y.levels <- sort(unique(y.value))
    if (is.logical(y.value) && !is.null(y.reference)) {
      y.reference <- as.logical(y.reference)
    }
  }
  q.total <- length(y.levels)
  if (q.total < 2L) {
    stop("Non-continuous families require at least two observed response levels.")
  }
  if (family == "binary" && q.total != 2L) {
    stop("The binary family requires exactly two observed response levels.")
  }
  if (family == "nominal") {
    if (is.null(y.reference)) {
      y.reference <- y.levels[1]
    }
    reference.index <- match(y.reference, y.levels)
    if (is.na(reference.index)) {
      stop("`y.reference` must match one of the observed response levels.")
    }
    q.set.index <- setdiff(seq_len(q.total), reference.index)
    q.set <- y.levels[q.set.index]
  } else if (family == "binary") {
    reference.index <- 1L
    y.reference <- y.levels[reference.index]
    q.set.index <- 2L
    q.set <- y.levels[q.set.index]
  } else {
    reference.index <- NA_integer_
    y.reference <- NULL
    q.set.index <- seq_len(q.total - 1L)
    q.set <- y.levels[q.set.index]
  }
  y.code <- match(y.value, y.levels)
  if (anyNA(y.code)) {
    stop("Unable to encode the response levels.")
  }
  y.by.q.vector <- lapply(seq_along(q.set.index), function(q) {
    if (family %in% c("binary", "nominal")) {
      as.integer(y.code == q.set.index[q])
    } else {
      as.integer(y.code <= q.set.index[q])
    }
  })
  list(
    y.input = y.value,
    y.levels = y.levels,
    q.total = q.total,
    n.q = length(q.set.index),
    q.set.index = q.set.index,
    q.set = q.set,
    y.reference = y.reference,
    y.by.q.vector = y.by.q.vector
  )
}
boostmtree.get.mu <- function(linear.predictor, family) {
  if (is.list(linear.predictor)) {
    return(lapply(linear.predictor, boostmtree.get.mu, family = family))
  }
  if (family == "continuous") {
    linear.predictor
  } else if (family %in% c("binary", "ordinal")) {
    plogis(linear.predictor)
  } else if (family == "nominal") {
    exp(linear.predictor)
  } else {
    stop("Unsupported family: ", family)
  }
}
boostmtree.get.mu.lambda <- function(linear.predictor, family) {
  if (is.list(linear.predictor)) {
    return(lapply(linear.predictor, boostmtree.get.mu.lambda, family = family))
  }
  if (family == "continuous") {
    rep(1, length(linear.predictor))
  } else if (family %in% c("binary", "ordinal")) {
    plogis(linear.predictor)
  } else if (family == "nominal") {
    exp(linear.predictor)
  } else {
    stop("Unsupported family: ", family)
  }
}
boostmtree.transform.h <- function(mu, family) {
  if (is.list(mu)) {
    return(lapply(mu, boostmtree.transform.h, family = family))
  }
  if (family == "continuous") {
    boostmtree.diag.matrix(rep(1, length(mu)))
  } else if (family %in% c("binary", "ordinal")) {
    boostmtree.diag.matrix(mu * (1 - mu))
  } else if (family == "nominal") {
    boostmtree.diag.matrix(mu)
  } else {
    stop("Unsupported family: ", family)
  }
}
boostmtree.build.prob.class <- function(mu, family, y.levels, q.set.index) {
  if (family == "continuous") {
    return(NULL)
  }
  n.subject <- length(mu[[1]])
  class.prob <- vector("list", length(y.levels))
  names(class.prob) <- as.character(y.levels)
  if (family == "binary") {
    class.prob[[1]] <- lapply(seq_len(n.subject), function(i) {
      pmax(0, pmin(1, 1 - mu[[1]][[i]]))
    })
    class.prob[[2]] <- mu[[1]]
    return(class.prob)
  }
  if (family == "nominal") {
    reference.prob <- lapply(seq_len(n.subject), function(i) {
      pmax(0, pmin(1, 1 - Reduce("+", lapply(seq_along(mu), function(q) mu[[q]][[i]]))))
    })
    non.reference.cursor <- 1L
    for (level.index in seq_along(y.levels)) {
      if (level.index %in% q.set.index) {
        class.prob[[level.index]] <- mu[[non.reference.cursor]]
        non.reference.cursor <- non.reference.cursor + 1L
      } else {
        class.prob[[level.index]] <- reference.prob
      }
    }
    return(class.prob)
  }
  class.prob[[1]] <- mu[[1]]
  if (length(mu) > 1L) {
    for (q in 2:length(mu)) {
      class.prob[[q]] <- lapply(seq_len(n.subject), function(i) {
        pmax(0, pmin(1, mu[[q]][[i]] - mu[[q - 1L]][[i]]))
      })
    }
  }
  class.prob[[length(y.levels)]] <- lapply(seq_len(n.subject), function(i) {
    pmax(0, pmin(1, 1 - mu[[length(mu)]][[i]]))
  })
  class.prob
}
