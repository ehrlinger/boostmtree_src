boostmtree <- function(
  x,
  tm = NULL,
  id = NULL,
  y,
  family = c("continuous", "binary", "nominal", "ordinal"),
  y.reference = NULL,
  M = 200,
  nu = 0.05,
  na.action = c("na.omit", "na.impute")[2],
  k = 5,
  mtry = NULL,
  n.knots = 10,
  d = 3,
  pen.ord = 3,
  lambda = NULL,
  rho = NULL,
  lambda.max = 1e6,
  lambda.iter = 2,
  svd.tol = 1e-6,
  verbose = TRUE,
  cv.flag = FALSE,
  eps = 1e-5,
  mod.grad = TRUE,
  nr.iter = 3,
  control = boostmtree.control()
) {
  family <- boostmtree.match.family(family)
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  control <- boostmtree.normalize.control(control)
  if (length(M) != 1L || is.na(M)) {
    stop("`M` must be a single positive integer.")
  }
  M <- as.integer(M)
  if (M < 1L) {
    stop("`M` must be at least 1.")
  }
  if (length(d) != 1L || is.na(d)) {
    stop("`d` must be a single integer.")
  }
  d <- as.integer(d)
  if (length(pen.ord) != 1L || is.na(pen.ord)) {
    stop("`pen.ord` must be a single non-negative integer.")
  }
  pen.ord <- as.integer(pen.ord)
  if (pen.ord < 0L) {
    stop("`pen.ord` must be non-negative.")
  }
  if (length(lambda.max) != 1L || is.na(lambda.max) || lambda.max <= 0) {
    stop("`lambda.max` must be a single positive numeric value.")
  }
  lambda.max <- as.numeric(lambda.max)
  if (length(svd.tol) != 1L || is.na(svd.tol) || svd.tol <= 0) {
    stop("`svd.tol` must be a single positive numeric value.")
  }
  svd.tol <- as.numeric(svd.tol)
  if (length(eps) != 1L || is.na(eps) || eps <= 0) {
    stop("`eps` must be a single positive numeric value.")
  }
  eps <- as.numeric(eps)
  verbose <- isTRUE(verbose)
  cv.flag <- isTRUE(cv.flag)
  mod.grad <- isTRUE(mod.grad)
  # ---------------------------------------------------------------------------
  # 1. Preprocess the subject-level and longitudinal inputs.
  # ---------------------------------------------------------------------------
  preprocessed <- boostmtree.preprocess.data(
    x = x,
    tm = tm,
    id = id,
    y = y,
    na.action = na.action
  )
  if (preprocessed$univariate) {
    mod.grad <- FALSE
    d <- -1L
  }
  # ---------------------------------------------------------------------------
  # 1b. Enforce an OOB-producing resampling scheme when OOB CV is requested.
  # ---------------------------------------------------------------------------
  control <- boostmtree.enforce.oob.control(
    control = control,
    cv.flag = cv.flag,
    n = preprocessed$n,
    M = M
  )
  # ---------------------------------------------------------------------------
  # 2. Encode the response in the canonical lower-case family vocabulary.
  # ---------------------------------------------------------------------------
  response.info <- boostmtree.prepare.response(
    y = preprocessed$y,
    family = family,
    y.reference = y.reference
  )
  # ---------------------------------------------------------------------------
  # 3. Build the time basis and the subject-specific design matrices.
  # ---------------------------------------------------------------------------
  time.design.info <- boostmtree.build.time.design(
    time.by.subject = preprocessed$time.by.subject,
    d = d,
    n.knots = n.knots
  )
  if (family == "continuous") {
    y.mean <- mean(preprocessed$y, na.rm = TRUE)
    y.sd <- sd(preprocessed$y, na.rm = TRUE)
    if (y.sd < 1e-6) {
      y.sd <- 1
    }
  } else {
    y.mean <- 0
    y.sd <- 1
  }
  y.org <- lapply(seq_len(response.info$n.q), function(q) {
    boostmtree.split.by.subject(
      response.info$y.by.q.vector[[q]],
      id = preprocessed$id,
      id.unique = preprocessed$id.unique
    )
  })
  # ---------------------------------------------------------------------------
  # 4. Fit the first-pass tree learner (`ntree = 1`) with cleaned controls.
  # ---------------------------------------------------------------------------
  fit.info <- boostmtree.fit.tree(list(
    x.subject = preprocessed$x.subject,
    x.long = preprocessed$x.long,
    x.var.names = preprocessed$x.var.names,
    time = preprocessed$time,
    time.by.subject = preprocessed$time.by.subject,
    time.unique = time.design.info$time.unique,
    id = preprocessed$id,
    id.unique = preprocessed$id.unique,
    y = preprocessed$y,
    y.by.subject = preprocessed$y.by.subject,
    y.by.q.vector = response.info$y.by.q.vector,
    y.org = y.org,
    y.mean = y.mean,
    y.sd = y.sd,
    y.levels = response.info$y.levels,
    y.reference = response.info$y.reference,
    family = family,
    n = preprocessed$n,
    ni = preprocessed$ni,
    n.q = response.info$n.q,
    q.total = response.info$q.total,
    q.set = response.info$q.set,
    q.set.index = response.info$q.set.index,
    x.tm = time.design.info$x.tm,
    time.design = time.design.info$time.design,
    time.unique = time.design.info$time.unique,
    df.time.design = time.design.info$df.time.design,
    d = time.design.info$d,
    pen.ord = pen.ord,
    M = M,
    nu = nu,
    k = k,
    mtry = mtry,
    lambda = lambda,
    rho = rho,
    lambda.max = lambda.max,
    lambda.iter = lambda.iter,
    svd.tol = svd.tol,
    verbose = verbose,
    cv.flag = cv.flag,
    eps = eps,
    mod.grad = mod.grad,
    nr.iter = nr.iter,
    control = control,
    na.action = na.action,
    univariate = preprocessed$univariate
  ))
  # ---------------------------------------------------------------------------
  # 5. Assemble the canonical return object.
  # ---------------------------------------------------------------------------
  object <- boostmtree.build.object(
    model.info = list(
      x.subject = preprocessed$x.subject,
      x.var.names = preprocessed$x.var.names,
      time.by.subject = preprocessed$time.by.subject,
      time.unique = time.design.info$time.unique,
      id = preprocessed$id,
      id.unique = preprocessed$id.unique,
      y.by.subject = preprocessed$y.by.subject,
      y.mean = y.mean,
      y.sd = y.sd,
      y.levels = response.info$y.levels,
      y.reference = response.info$y.reference,
      family = family,
      na.action = na.action,
      n = preprocessed$n,
      ni = preprocessed$ni,
      n.q = response.info$n.q,
      q.total = response.info$q.total,
      q.set = response.info$q.set,
      q.set.index = response.info$q.set.index,
      x.tm = time.design.info$x.tm,
      time.design = time.design.info$time.design,
      df.time.design = time.design.info$df.time.design,
      d = time.design.info$d,
      pen.ord = pen.ord,
      k = k,
      M = M,
      nu = if (length(nu) == 1L) rep(nu, 2L) else nu,
      control = control,
      cv.flag = cv.flag,
      univariate = preprocessed$univariate
    ),
    fit.info = fit.info
  )
  invisible(object)
}
