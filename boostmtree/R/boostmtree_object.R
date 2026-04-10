boostmtree.flatten.single.response <- function(x, family) {
  if (is.null(x)) {
    return(NULL)
  }
  if (family %in% c("nominal", "ordinal")) {
    return(x)
  }
  if (is.list(x) && length(x) == 1L) {
    return(x[[1L]])
  }
  x
}
boostmtree.flatten.single.path <- function(x, family) {
  if (is.null(x)) {
    return(NULL)
  }
  if (family %in% c("nominal", "ordinal")) {
    return(x)
  }
  as.vector(x)
}
boostmtree.flatten.oob.count <- function(x, n.q) {
  if (is.null(x)) {
    return(NULL)
  }
  if (n.q == 1L) {
    return(as.integer(x[, 1L]))
  }
  x
}
boostmtree.flatten.single.error <- function(x, family, y.sd) {
  if (is.null(x)) {
    return(NULL)
  }
  if (family %in% c("nominal", "ordinal")) {
    return(lapply(seq_along(x), function(q) x[[q]] / y.sd))
  }
  x[[1L]] / y.sd
}
boostmtree.build.object <- function(model.info, fit.info) {
  prob.class <- boostmtree.build.prob.class(
    mu = fit.info$mu,
    family = model.info$family,
    y.levels = model.info$y.levels,
    q.set.index = model.info$q.set.index
  )
  learner.used <- if (model.info$df.time.design == 1L) {
    "tree.learner"
  } else {
    "mtree.pspline.learner"
  }
  object <- list(
    x = model.info$x.subject,
    x.var.names = model.info$x.var.names,
    time = model.info$time.by.subject,
    time.unique = model.info$time.unique,
    id = model.info$id,
    id.unique = model.info$id.unique,
    y = model.info$y.by.subject,
    y.org = boostmtree.flatten.single.response(fit.info$y.org, family = model.info$family),
    family = model.info$family,
    y.mean = model.info$y.mean,
    y.sd = model.info$y.sd,
    y.levels = model.info$y.levels,
    y.reference = model.info$y.reference,
    na.action = model.info$na.action,
    n = model.info$n,
    ni = model.info$ni,
    n.q = model.info$n.q,
    q.total = model.info$q.total,
    q.set = model.info$q.set,
    mu = boostmtree.flatten.single.response(fit.info$mu, family = model.info$family),
    prob.class = prob.class,
    lambda = boostmtree.flatten.single.path(fit.info$lambda, family = model.info$family),
    phi = boostmtree.flatten.single.path(fit.info$phi, family = model.info$family),
    rho = boostmtree.flatten.single.path(fit.info$rho, family = model.info$family),
    gamma = fit.info$gamma,
    base.learner = fit.info$base.learner,
    membership = fit.info$membership,
    x.tm = model.info$x.tm,
    time.design = model.info$time.design,
    d = model.info$d,
    pen.ord = model.info$pen.ord,
    k = model.info$k,
    M = model.info$M,
    nu = model.info$nu,
    ntree = model.info$control$ntree,
    control = model.info$control,
    cv.flag = model.info$cv.flag,
    err.rate = boostmtree.flatten.single.error(
      fit.info$err.rate,
      family = model.info$family,
      y.sd = model.info$y.sd
    ),
    rmse = if (!is.null(fit.info$rmse)) {
      as.numeric(fit.info$rmse / model.info$y.sd)
    } else {
      NULL
    },
    m.opt = fit.info$m.opt,
    gamma.i.list = fit.info$gamma.i.list,
    oob.available = isTRUE(fit.info$oob.available),
    oob.subject.count = boostmtree.flatten.oob.count(
      fit.info$oob.subject.count,
      n.q = model.info$n.q
    ),
    univariate = model.info$univariate
  )
  class(object) <- c("boostmtree", "grow", learner.used)
  invisible(object)
}
