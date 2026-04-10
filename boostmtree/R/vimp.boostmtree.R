boostmtree.vimp.is.grow.object <- function(object) {
  all(c("boostmtree", "grow") %in% class(object))
}
boostmtree.vimp.is.predict.object <- function(object) {
  all(c("boostmtree", "predict") %in% class(object))
}
boostmtree.vimp.field <- function(x, canonical, legacy = NULL, default = NULL) {
  if (!is.null(x[[canonical]])) {
    return(x[[canonical]])
  }
  if (!is.null(legacy) && !is.null(x[[legacy]])) {
    return(x[[legacy]])
  }
  default
}
boostmtree.vimp.split.by.subject <- function(values, id, id.unique) {
  unname(split(values, factor(id, levels = id.unique)))
}
boostmtree.vimp.as.q.list <- function(x, n.q) {
  if (is.null(x)) {
    return(NULL)
  }
  if (n.q == 1L) {
    if (is.matrix(x)) {
      return(list(x))
    }
    if (!is.list(x)) {
      return(list(x))
    }
    if (length(x) == 0L) {
      return(list(x))
    }
    if (!is.list(x[[1L]])) {
      return(list(x))
    }
  }
  x
}
boostmtree.vimp.select.variables <- function(x.var.names, x.names = NULL, joint = FALSE) {
  p.total <- length(x.var.names)
  if (is.null(x.names)) {
    selected.index <- seq_len(p.total)
    selected.names <- x.var.names
  } else {
    selected.index <- match(x.names, x.var.names)
    if (anyNA(selected.index)) {
      stop("`x.names` must match names from the fitted object.")
    }
    if (joint) {
      selected.names <- "joint.vimp"
    } else {
      selected.names <- x.names
    }
  }
  list(
    selected.index = selected.index,
    selected.names = selected.names,
    p = if (joint) 1L else length(selected.index)
  )
}
boostmtree.vimp.get.time.design <- function(object) {
  time.design <- boostmtree.vimp.field(object, "time.design", legacy = "D")
  if (is.null(time.design)) {
    stop("Unable to find the time-design matrices in the supplied object.")
  }
  time.design
}
boostmtree.vimp.get.df.time.design <- function(object) {
  x.tm <- boostmtree.vimp.field(object, "x.tm", legacy = "X.tm")
  if (!is.null(x.tm)) {
    return(ncol(x.tm))
  }
  ncol(boostmtree.vimp.get.time.design(object)[[1L]])
}
boostmtree.vimp.get.nu.vec <- function(object, df.time.design) {
  nu.vec <- boostmtree.vimp.field(object, "nu.vec")
  if (!is.null(nu.vec)) {
    return(as.numeric(nu.vec))
  }
  nu <- boostmtree.vimp.field(object, "nu")
  if (is.null(nu)) {
    stop("Unable to recover the learning-rate vector from the supplied object.")
  }
  if (length(nu) == 1L) {
    rep(as.numeric(nu), df.time.design)
  } else {
    c(as.numeric(nu[1L]), rep(as.numeric(nu[2L]), df.time.design - 1L))
  }
}
boostmtree.vimp.q.names <- function(q.set, n.q, family) {
  if (n.q == 1L && identical(family, "continuous")) {
    return("response")
  }
  if (is.null(q.set) || all(is.na(q.set))) {
    paste0("q", seq_len(n.q))
  } else {
    as.character(q.set)
  }
}
boostmtree.vimp.parse.gamma.step <- function(gamma.step) {
  if (is.null(gamma.step)) {
    stop("Encountered a NULL gamma step.")
  }
  if (is.list(gamma.step) && !is.null(gamma.step$node.label)) {
    coefficients <- gamma.step$coefficients
    if (is.null(dim(coefficients))) {
      coefficients <- matrix(coefficients, nrow = 1L)
    }
    return(list(node.label = gamma.step$node.label, coefficients = coefficients))
  }
  if (is.matrix(gamma.step)) {
    return(list(node.label = gamma.step[, 1L, drop = TRUE], coefficients = gamma.step[, -1L, drop = FALSE]))
  }
  stop("Unsupported gamma representation.")
}
boostmtree.vimp.extract.coefficients <- function(gamma.step, node.label) {
  parsed.gamma <- boostmtree.vimp.parse.gamma.step(gamma.step)
  row.index <- match(as.character(node.label), as.character(parsed.gamma$node.label))
  if (is.na(row.index)) {
    stop("Unable to match a node label in the gamma step.")
  }
  as.numeric(parsed.gamma$coefficients[row.index, , drop = FALSE])
}
boostmtree.vimp.relative.increase <- function(error.value, baseline.error) {
  if (is.na(baseline.error) || abs(baseline.error) < 1e-12) {
    return(NA_real_)
  }
  (error.value - baseline.error) / baseline.error
}
boostmtree.vimp.oob.count <- function(object, m.opt, n.q) {
  oob.subject.count <- boostmtree.vimp.field(object, "oob.subject.count")
  if (!is.null(oob.subject.count)) {
    if (n.q == 1L && is.null(dim(oob.subject.count))) {
      return(matrix(as.integer(oob.subject.count), ncol = 1L))
    }
    return(as.matrix(oob.subject.count))
  }
  base.learner <- boostmtree.vimp.field(object, "base.learner", legacy = "baselearner")
  if (is.null(base.learner)) {
    return(NULL)
  }
  count.matrix <- matrix(NA_integer_, nrow = max(m.opt), ncol = n.q)
  for (q in seq_len(n.q)) {
    for (m in seq_len(m.opt[q])) {
      count.matrix[m, q] <- sum(base.learner[[q]][[m]]$inbag == 0L)
    }
  }
  count.matrix
}
boostmtree.vimp.require.oob <- function(object, m.opt, n.q) {
  oob.count <- boostmtree.vimp.oob.count(object = object, m.opt = m.opt, n.q = n.q)
  if (is.null(oob.count)) {
    return(invisible(TRUE))
  }
  bad.iteration <- lapply(seq_len(n.q), function(q) {
    which(oob.count[seq_len(m.opt[q]), q] <= 0L)
  })
  if (all(vapply(bad.iteration, length, integer(1)) == 0L)) {
    return(invisible(TRUE))
  }
  bad.message <- paste(vapply(seq_len(n.q), function(q) {
    if (length(bad.iteration[[q]]) == 0L) {
      return(character(1))
    }
    paste0("q", q, ": ", paste(bad.iteration[[q]], collapse = ", "))
  }, character(1)), collapse = "; ")
  bad.message <- gsub("(^; |; $|^$)", "", bad.message)
  stop(
    "Grow-object variable importance requires out-of-bag subjects for every boosting iteration used in the importance calculation. ",
    "Refit with `cv.flag = TRUE` and an OOB-producing resampling scheme, such as the default `control = boostmtree.control()` or `control = boostmtree.control(seed = ...)`. ",
    if (nzchar(bad.message)) paste0("Iterations with no OOB subjects: ", bad.message, ".") else ""
  )
}
boostmtree.vimp.response.by.q <- function(object, model, family, n.q) {
  response.by.q <- boostmtree.vimp.field(object, "y.org", legacy = "Yorg")
  if (!is.null(response.by.q)) {
    return(boostmtree.vimp.as.q.list(response.by.q, n.q))
  }
  response.by.q <- boostmtree.vimp.field(object, "Y")
  if (!is.null(response.by.q)) {
    return(boostmtree.vimp.as.q.list(response.by.q, n.q))
  }
  raw.y <- boostmtree.vimp.field(object, "y")
  if (is.null(raw.y)) {
    return(NULL)
  }
  if (is.list(raw.y)) {
    if (n.q == 1L) {
      return(list(raw.y))
    }
    return(raw.y)
  }
  id <- boostmtree.vimp.field(object, "id", default = boostmtree.vimp.field(model, "id"))
  id.unique <- boostmtree.vimp.field(object, "id.unique", default = boostmtree.vimp.field(model, "id.unique", default = sort(unique(id))))
  if (family == "continuous") {
    return(list(boostmtree.vimp.split.by.subject(raw.y, id, id.unique)))
  }
  y.levels <- boostmtree.vimp.field(model, "y.levels", legacy = "y.unq")
  q.set <- boostmtree.vimp.field(model, "q.set", legacy = "Q_set")
  if (is.factor(raw.y)) {
    raw.value <- as.character(droplevels(raw.y))
  } else {
    raw.value <- raw.y
  }
  response.index <- match(as.character(raw.value), as.character(y.levels))
  if (anyNA(response.index)) {
    stop("Response levels in the prediction object do not match the fitted model.")
  }
  q.set.index <- match(as.character(q.set), as.character(y.levels))
  y.by.q.vector <- lapply(seq_along(q.set.index), function(q) {
    if (family %in% c("binary", "nominal")) {
      as.integer(response.index == q.set.index[q])
    } else {
      as.integer(response.index <= q.set.index[q])
    }
  })
  lapply(y.by.q.vector, function(z) boostmtree.vimp.split.by.subject(z, id, id.unique))
}
boostmtree.vimp.predict.membership <- function(base.learner, newdata, k, na.action = "na.impute") {
  c(predict.rfsrc(base.learner, newdata = newdata, membership = TRUE, ptn.count = k, importance = "none", na.action = na.action)$ptn.membership)
}
boostmtree.vimp.make.reference <- function(l.pred.db.vimp, family, n.q, ni, n, p) {
  l.pred.ref.vimp <- vector("list", p)
  for (k.index in seq_len(p)) {
    if (family == "nominal") {
      l.pred.ref.vimp[[k.index]] <- lapply(seq_len(n), function(i) {
        ref.main <- log((1 + Reduce("+", lapply(seq_len(n.q), function(q) exp(l.pred.db.vimp[[q]][[k.index]][[i]]$main))))^(-1))
        ref.interaction <- if (!is.null(l.pred.db.vimp[[1L]][[k.index]][[i]]$interaction)) {
          log((1 + Reduce("+", lapply(seq_len(n.q), function(q) exp(l.pred.db.vimp[[q]][[k.index]][[i]]$interaction))))^(-1))
        } else {
          NULL
        }
        ref.time <- if (!is.null(l.pred.db.vimp[[1L]][[k.index]][[i]]$time)) {
          log((1 + Reduce("+", lapply(seq_len(n.q), function(q) exp(l.pred.db.vimp[[q]][[k.index]][[i]]$time))))^(-1))
        } else {
          NULL
        }
        list(main = ref.main, interaction = ref.interaction, time = ref.time)
      })
    } else {
      l.pred.ref.vimp[[k.index]] <- lapply(seq_len(n), function(i) {
        list(
          main = rep(0, ni[i]),
          interaction = if (!is.null(l.pred.db.vimp[[1L]][[k.index]][[i]]$interaction)) rep(0, ni[i]) else NULL,
          time = if (!is.null(l.pred.db.vimp[[1L]][[k.index]][[i]]$time)) rep(0, ni[i]) else NULL
        )
      })
    }
  }
  l.pred.ref.vimp
}
boostmtree.vimp.build.object <- function(main, interaction, time.effect, x.var.names, q.set, family, source, joint, metric, baseline.rmse, m.opt) {
  n.q <- ncol(main)
  q.names <- boostmtree.vimp.q.names(q.set = q.set, n.q = n.q, family = family)
  colnames(main) <- q.names
  if (!is.null(interaction)) {
    colnames(interaction) <- q.names
  }
  if (!is.null(time.effect)) {
    names(time.effect) <- q.names
  }
  object <- list(
    main = main,
    interaction = interaction,
    time.effect = time.effect,
    x.var.names = x.var.names,
    q.set = q.set,
    family = family,
    source = source,
    joint = joint,
    metric = metric,
    baseline.rmse = baseline.rmse,
    m.opt = m.opt
  )
  class(object) <- c("vimp.boostmtree", "boostmtree.vimp")
  invisible(object)
}
boostmtree.vimp.from.grow <- function(object, x.names = NULL, joint = FALSE) {
  if (!isTRUE(boostmtree.vimp.field(object, "cv.flag"))) {
    stop("Grow-object variable importance requires `cv.flag = TRUE`.")
  }
  x <- boostmtree.vimp.field(object, "x")
  x.var.names <- boostmtree.vimp.field(object, "x.var.names", legacy = "xvar.names")
  variable.info <- boostmtree.vimp.select.variables(x.var.names = x.var.names, x.names = x.names, joint = joint)
  y.mean <- boostmtree.vimp.field(object, "y.mean", legacy = "ymean")
  y.sd <- boostmtree.vimp.field(object, "y.sd", legacy = "ysd")
  df.time.design <- boostmtree.vimp.get.df.time.design(object)
  m.opt <- as.integer(boostmtree.vimp.field(object, "m.opt", legacy = "Mopt"))
  ni <- boostmtree.vimp.field(object, "ni")
  k <- as.integer(boostmtree.vimp.field(object, "k", legacy = "K"))
  n.q <- as.integer(boostmtree.vimp.field(object, "n.q", legacy = "n.Q"))
  q.set <- boostmtree.vimp.field(object, "q.set", legacy = "Q_set")
  gamma.i.list <- boostmtree.vimp.field(object, "gamma.i.list")
  membership <- boostmtree.vimp.field(object, "membership")
  time.design <- boostmtree.vimp.get.time.design(object)
  n <- as.integer(boostmtree.vimp.field(object, "n"))
  rmse <- boostmtree.vimp.field(object, "rmse") * y.sd
  nu.vec <- boostmtree.vimp.get.nu.vec(object, df.time.design)
  family <- boostmtree.vimp.field(object, "family")
  y.by.q <- boostmtree.vimp.response.by.q(object, object, family = family, n.q = n.q)
  base.learner <- boostmtree.vimp.field(object, "base.learner", legacy = "baselearner")
  if (is.null(gamma.i.list)) {
    stop("Grow-object variable importance requires `gamma.i.list` from the CV fit path.")
  }
  if (anyNA(m.opt)) {
    stop("Unable to determine the optimized number of iterations from the fitted object.")
  }
  boostmtree.vimp.require.oob(object = object, m.opt = m.opt, n.q = n.q)
  oob.list <- lapply(seq_len(n.q), function(q) vector("list", m.opt[q]))
  membership.noise.list <- vector("list", n.q)
  for (q in seq_len(n.q)) {
    membership.noise.list[[q]] <- lapply(seq_len(m.opt[q]), function(m) {
      oob <- which(base.learner[[q]][[m]]$inbag == 0)
      oob.list[[q]][[m]] <<- oob
      n.oob <- length(oob)
      if (n.oob == 0L) {
        return(matrix(NA, nrow = 0L, ncol = variable.info$p))
      }
      x.noise <- do.call(rbind, lapply(seq_len(variable.info$p), function(k.index) {
        x.k <- x[oob, , drop = FALSE]
        if (joint) {
          x.k[, variable.info$selected.index] <- x.k[sample.int(nrow(x.k)), variable.info$selected.index, drop = FALSE]
        } else {
          x.k[, variable.info$selected.index[k.index]] <- sample(x.k[, variable.info$selected.index[k.index]], size = nrow(x.k), replace = FALSE)
        }
        x.k
      }))
      membership.noise <- boostmtree.vimp.predict.membership(base.learner = base.learner[[q]][[m]], newdata = x.noise, k = k, na.action = "na.impute")
      matrix(membership.noise, nrow = n.oob, byrow = FALSE)
    })
  }
  l.pred.db.vimp <- lapply(seq_len(n.q), function(q) vector("list", variable.info$p))
  if (df.time.design > 1L) {
    for (q in seq_len(n.q)) {
      for (k.index in seq_len(variable.info$p)) {
        l.pred.db.vimp[[q]][[k.index]] <- lapply(seq_len(n), function(i) {
          l.pred.db.main.i <- rep(0, ni[i])
          l.pred.db.interaction.i <- rep(0, ni[i])
          l.pred.db.time.i <- if (k.index == variable.info$p) rep(0, ni[i]) else NULL
          for (m in seq_len(m.opt[q])) {
            if (i %in% oob.list[[q]][[m]]) {
              noise.label <- membership.noise.list[[q]][[m]][which(oob.list[[q]][[m]] == i), k.index]
              gamma.step.i <- gamma.i.list[[q]][[m]][[i]]
              org.label <- membership[[q]][[m]][i]
              gamma.noise.i <- boostmtree.vimp.extract.coefficients(gamma.step.i, noise.label)
              gamma.org.i <- boostmtree.vimp.extract.coefficients(gamma.step.i, org.label)
              gamma.main <- c(gamma.noise.i[1L], gamma.org.i[-1L])
              gamma.interaction <- c(gamma.org.i[1L], gamma.noise.i[-1L])
              l.pred.db.main.i <- l.pred.db.main.i + c(time.design[[i]] %*% (gamma.main * nu.vec))
              l.pred.db.interaction.i <- l.pred.db.interaction.i + c(time.design[[i]] %*% (gamma.interaction * nu.vec))
              if (k.index == variable.info$p) {
                n.d <- nrow(time.design[[i]])
                l.pred.db.time.i <- l.pred.db.time.i + c(time.design[[i]][sample.int(n.d, n.d, replace = TRUE), , drop = FALSE] %*% (gamma.org.i * nu.vec))
              }
            }
          }
          list(main = l.pred.db.main.i, interaction = l.pred.db.interaction.i, time = l.pred.db.time.i)
        })
      }
    }
    l.pred.ref.vimp <- boostmtree.vimp.make.reference(l.pred.db.vimp = l.pred.db.vimp, family = family, n.q = n.q, ni = ni, n = n, p = variable.info$p)
    vimp.main <- vimp.interaction <- matrix(NA_real_, nrow = variable.info$p, ncol = n.q)
    vimp.time <- rep(NA_real_, n.q)
    for (q in seq_len(n.q)) {
      for (k.index in seq_len(variable.info$p)) {
        mu.main <- lapply(seq_len(n), function(i) boostmtree.get.mu(linear.predictor = (l.pred.db.vimp[[q]][[k.index]][[i]]$main + l.pred.ref.vimp[[k.index]][[i]]$main) * y.sd + y.mean, family = family))
        mu.interaction <- lapply(seq_len(n), function(i) boostmtree.get.mu(linear.predictor = (l.pred.db.vimp[[q]][[k.index]][[i]]$interaction + l.pred.ref.vimp[[k.index]][[i]]$interaction) * y.sd + y.mean, family = family))
        err.main <- boostmtree.l2.dist(y.by.q[[q]], mu.main)
        err.interaction <- boostmtree.l2.dist(y.by.q[[q]], mu.interaction)
        vimp.main[k.index, q] <- boostmtree.vimp.relative.increase(err.main, rmse[q])
        vimp.interaction[k.index, q] <- boostmtree.vimp.relative.increase(err.interaction, rmse[q])
        if (k.index == variable.info$p) {
          mu.time <- lapply(seq_len(n), function(i) boostmtree.get.mu(linear.predictor = (l.pred.db.vimp[[q]][[k.index]][[i]]$time + l.pred.ref.vimp[[k.index]][[i]]$time) * y.sd + y.mean, family = family))
          err.time <- boostmtree.l2.dist(y.by.q[[q]], mu.time)
          vimp.time[q] <- boostmtree.vimp.relative.increase(err.time, rmse[q])
        }
      }
    }
    rownames(vimp.main) <- variable.info$selected.names
    rownames(vimp.interaction) <- paste(variable.info$selected.names, "time", sep = ":")
    return(boostmtree.vimp.build.object(main = vimp.main, interaction = vimp.interaction, time.effect = vimp.time, x.var.names = variable.info$selected.names, q.set = q.set, family = family, source = "grow", joint = joint, metric = "relative increase in OOB RMSE", baseline.rmse = rmse, m.opt = m.opt))
  }
  for (q in seq_len(n.q)) {
    for (k.index in seq_len(variable.info$p)) {
      l.pred.db.vimp[[q]][[k.index]] <- lapply(seq_len(n), function(i) {
        l.pred.db.i <- rep(0, ni[i])
        for (m in seq_len(m.opt[q])) {
          if (i %in% oob.list[[q]][[m]]) {
            noise.label <- membership.noise.list[[q]][[m]][which(oob.list[[q]][[m]] == i), k.index]
            gamma.step.i <- gamma.i.list[[q]][[m]][[i]]
            gamma.noise.i <- boostmtree.vimp.extract.coefficients(gamma.step.i, noise.label)
            l.pred.db.i <- l.pred.db.i + c(time.design[[i]] %*% (gamma.noise.i * nu.vec))
          }
        }
        list(main = l.pred.db.i, interaction = NULL, time = NULL)
      })
    }
  }
  l.pred.ref.vimp <- boostmtree.vimp.make.reference(l.pred.db.vimp = l.pred.db.vimp, family = family, n.q = n.q, ni = ni, n = n, p = variable.info$p)
  vimp.main <- matrix(NA_real_, nrow = variable.info$p, ncol = n.q)
  for (q in seq_len(n.q)) {
    for (k.index in seq_len(variable.info$p)) {
      mu.main <- lapply(seq_len(n), function(i) boostmtree.get.mu(linear.predictor = (l.pred.db.vimp[[q]][[k.index]][[i]]$main + l.pred.ref.vimp[[k.index]][[i]]$main) * y.sd + y.mean, family = family))
      err.main <- boostmtree.l2.dist(y.by.q[[q]], mu.main)
      vimp.main[k.index, q] <- boostmtree.vimp.relative.increase(err.main, rmse[q])
    }
  }
  rownames(vimp.main) <- variable.info$selected.names
  boostmtree.vimp.build.object(main = vimp.main, interaction = NULL, time.effect = NULL, x.var.names = variable.info$selected.names, q.set = q.set, family = family, source = "grow", joint = joint, metric = "relative increase in OOB RMSE", baseline.rmse = rmse, m.opt = m.opt)
}
boostmtree.vimp.from.predict <- function(object, x.names = NULL, joint = FALSE) {
  model <- boostmtree.vimp.field(object, "boost.obj")
  if (is.null(model)) {
    stop("Predict objects must contain a `boost.obj` component.")
  }
  x <- boostmtree.vimp.field(object, "x")
  x.var.names <- boostmtree.vimp.field(object, "x.var.names", legacy = "xvar.names")
  variable.info <- boostmtree.vimp.select.variables(x.var.names = x.var.names, x.names = x.names, joint = joint)
  if (is.null(boostmtree.vimp.field(object, "y")) && is.null(boostmtree.vimp.field(object, "Y")) && is.null(boostmtree.vimp.field(object, "y.org", legacy = "Yorg"))) {
    stop("Prediction-object variable importance requires observed responses in the prediction object.")
  }
  family <- boostmtree.vimp.field(object, "family", default = boostmtree.vimp.field(model, "family"))
  y.mean <- boostmtree.vimp.field(object, "y.mean", legacy = "ymean", default = boostmtree.vimp.field(model, "y.mean", legacy = "ymean"))
  y.sd <- boostmtree.vimp.field(object, "y.sd", legacy = "ysd", default = boostmtree.vimp.field(model, "y.sd", legacy = "ysd"))
  n <- as.integer(boostmtree.vimp.field(object, "n"))
  ni <- boostmtree.vimp.field(object, "ni")
  k <- as.integer(boostmtree.vimp.field(object, "k", legacy = "K", default = boostmtree.vimp.field(model, "k", legacy = "K")))
  n.q <- as.integer(boostmtree.vimp.field(object, "n.q", legacy = "n.Q", default = boostmtree.vimp.field(model, "n.q", legacy = "n.Q")))
  q.set <- boostmtree.vimp.field(object, "q.set", legacy = "Q_set", default = boostmtree.vimp.field(model, "q.set", legacy = "Q_set"))
  df.time.design <- boostmtree.vimp.get.df.time.design(model)
  time.design <- boostmtree.vimp.get.time.design(object)
  nu.vec <- boostmtree.vimp.get.nu.vec(object, df.time.design)
  m.opt <- as.integer(boostmtree.vimp.field(object, "m.opt", legacy = "Mopt"))
  base.learner <- boostmtree.vimp.field(object, "base.learner", legacy = "baselearner")
  gamma <- boostmtree.vimp.field(object, "gamma")
  membership <- boostmtree.vimp.field(object, "membership")
  rmse <- boostmtree.vimp.field(object, "rmse") * y.sd
  y.by.q <- boostmtree.vimp.response.by.q(object, model, family = family, n.q = n.q)
  if (anyNA(m.opt)) {
    stop("Unable to determine the optimized number of iterations for the prediction object.")
  }
  membership.noise <- vector("list", n.q)
  for (q in seq_len(n.q)) {
    membership.noise[[q]] <- lapply(seq_len(m.opt[q]), function(m) {
      x.noise <- do.call(rbind, lapply(seq_len(variable.info$p), function(k.index) {
        x.k <- x
        if (joint) {
          x.k[, variable.info$selected.index] <- x.k[sample.int(nrow(x.k)), variable.info$selected.index, drop = FALSE]
        } else {
          x.k[, variable.info$selected.index[k.index]] <- sample(x.k[, variable.info$selected.index[k.index]], size = nrow(x.k), replace = FALSE)
        }
        x.k
      }))
      boostmtree.vimp.predict.membership(base.learner = base.learner[[q]][[m]], newdata = x.noise, k = k, na.action = "na.impute")
    })
  }
  l.pred.db.vimp <- lapply(seq_len(n.q), function(q) vector("list", variable.info$p))
  if (df.time.design > 1L) {
    for (q in seq_len(n.q)) {
      for (k.index in seq_len(variable.info$p)) {
        l.pred.db.vimp[[q]][[k.index]] <- lapply(seq_len(n), function(i) {
          l.pred.db.main.i <- rep(0, ni[i])
          l.pred.db.interaction.i <- rep(0, ni[i])
          l.pred.db.time.i <- if (k.index == variable.info$p) rep(0, ni[i]) else NULL
          for (m in seq_len(m.opt[q])) {
            gamma.step <- gamma[[q]][[m]]
            org.label <- membership[[q]][[m]][i]
            membership.k <- membership.noise[[q]][[m]][((k.index - 1L) * n + 1L):(k.index * n)]
            noise.label <- membership.k[i]
            gamma.org.i <- boostmtree.vimp.extract.coefficients(gamma.step, org.label)
            gamma.noise.i <- boostmtree.vimp.extract.coefficients(gamma.step, noise.label)
            gamma.main <- c(gamma.noise.i[1L], gamma.org.i[-1L])
            gamma.interaction <- c(gamma.org.i[1L], gamma.noise.i[-1L])
            l.pred.db.main.i <- l.pred.db.main.i + c(time.design[[i]] %*% (gamma.main * nu.vec))
            l.pred.db.interaction.i <- l.pred.db.interaction.i + c(time.design[[i]] %*% (gamma.interaction * nu.vec))
            if (k.index == variable.info$p) {
              n.d <- nrow(time.design[[i]])
              l.pred.db.time.i <- l.pred.db.time.i + c(time.design[[i]][sample.int(n.d, n.d, replace = TRUE), , drop = FALSE] %*% (gamma.org.i * nu.vec))
            }
          }
          list(main = l.pred.db.main.i, interaction = l.pred.db.interaction.i, time = l.pred.db.time.i)
        })
      }
    }
    l.pred.ref.vimp <- boostmtree.vimp.make.reference(l.pred.db.vimp = l.pred.db.vimp, family = family, n.q = n.q, ni = ni, n = n, p = variable.info$p)
    vimp.main <- vimp.interaction <- matrix(NA_real_, nrow = variable.info$p, ncol = n.q)
    vimp.time <- rep(NA_real_, n.q)
    for (q in seq_len(n.q)) {
      for (k.index in seq_len(variable.info$p)) {
        mu.main <- lapply(seq_len(n), function(i) boostmtree.get.mu(linear.predictor = (l.pred.db.vimp[[q]][[k.index]][[i]]$main + l.pred.ref.vimp[[k.index]][[i]]$main) * y.sd + y.mean, family = family))
        mu.interaction <- lapply(seq_len(n), function(i) boostmtree.get.mu(linear.predictor = (l.pred.db.vimp[[q]][[k.index]][[i]]$interaction + l.pred.ref.vimp[[k.index]][[i]]$interaction) * y.sd + y.mean, family = family))
        err.main <- boostmtree.l2.dist(y.by.q[[q]], mu.main)
        err.interaction <- boostmtree.l2.dist(y.by.q[[q]], mu.interaction)
        vimp.main[k.index, q] <- boostmtree.vimp.relative.increase(err.main, rmse[q])
        vimp.interaction[k.index, q] <- boostmtree.vimp.relative.increase(err.interaction, rmse[q])
        if (k.index == variable.info$p) {
          mu.time <- lapply(seq_len(n), function(i) boostmtree.get.mu(linear.predictor = (l.pred.db.vimp[[q]][[k.index]][[i]]$time + l.pred.ref.vimp[[k.index]][[i]]$time) * y.sd + y.mean, family = family))
          err.time <- boostmtree.l2.dist(y.by.q[[q]], mu.time)
          vimp.time[q] <- boostmtree.vimp.relative.increase(err.time, rmse[q])
        }
      }
    }
    rownames(vimp.main) <- variable.info$selected.names
    rownames(vimp.interaction) <- paste(variable.info$selected.names, "time", sep = ":")
    return(boostmtree.vimp.build.object(main = vimp.main, interaction = vimp.interaction, time.effect = vimp.time, x.var.names = variable.info$selected.names, q.set = q.set, family = family, source = "predict", joint = joint, metric = "relative increase in test-set RMSE", baseline.rmse = rmse, m.opt = m.opt))
  }
  for (q in seq_len(n.q)) {
    for (k.index in seq_len(variable.info$p)) {
      l.pred.db.vimp[[q]][[k.index]] <- lapply(seq_len(n), function(i) {
        l.pred.db.i <- rep(0, ni[i])
        for (m in seq_len(m.opt[q])) {
          gamma.step <- gamma[[q]][[m]]
          membership.k <- membership.noise[[q]][[m]][((k.index - 1L) * n + 1L):(k.index * n)]
          noise.label <- membership.k[i]
          gamma.noise.i <- boostmtree.vimp.extract.coefficients(gamma.step, noise.label)
          l.pred.db.i <- l.pred.db.i + c(time.design[[i]] %*% (gamma.noise.i * nu.vec))
        }
        list(main = l.pred.db.i, interaction = NULL, time = NULL)
      })
    }
  }
  l.pred.ref.vimp <- boostmtree.vimp.make.reference(l.pred.db.vimp = l.pred.db.vimp, family = family, n.q = n.q, ni = ni, n = n, p = variable.info$p)
  vimp.main <- matrix(NA_real_, nrow = variable.info$p, ncol = n.q)
  for (q in seq_len(n.q)) {
    for (k.index in seq_len(variable.info$p)) {
      mu.main <- lapply(seq_len(n), function(i) boostmtree.get.mu(linear.predictor = (l.pred.db.vimp[[q]][[k.index]][[i]]$main + l.pred.ref.vimp[[k.index]][[i]]$main) * y.sd + y.mean, family = family))
      err.main <- boostmtree.l2.dist(y.by.q[[q]], mu.main)
      vimp.main[k.index, q] <- boostmtree.vimp.relative.increase(err.main, rmse[q])
    }
  }
  rownames(vimp.main) <- variable.info$selected.names
  boostmtree.vimp.build.object(main = vimp.main, interaction = NULL, time.effect = NULL, x.var.names = variable.info$selected.names, q.set = q.set, family = family, source = "predict", joint = joint, metric = "relative increase in test-set RMSE", baseline.rmse = rmse, m.opt = m.opt)
}
vimp.boostmtree <- function(object, x.names = NULL, joint = FALSE) {
  if (!boostmtree.vimp.is.grow.object(object) && !boostmtree.vimp.is.predict.object(object)) {
    stop("This function only works for `(boostmtree, grow)` or `(boostmtree, predict)` objects.")
  }
  if (boostmtree.vimp.is.grow.object(object)) {
    return(boostmtree.vimp.from.grow(object = object, x.names = x.names, joint = joint))
  }
  boostmtree.vimp.from.predict(object = object, x.names = x.names, joint = joint)
}
