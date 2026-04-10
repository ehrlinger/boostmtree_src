# Keep only utility functions that are still plausibly needed by older, not-yet-
# refactored code. Remove the duplicated modular code and legacy hidden-option
# parsers.
.bt_match_family <- function(family) {
  family <- tolower(as.character(family)[1L])
  match.arg(family, c("continuous", "binary", "nominal", "ordinal"))
}
# Create a diagonal matrix from a vector.
DiagMat <- function(x) {
  if (length(x) == 1L) {
    matrix(x, nrow = 1L, ncol = 1L)
  } else {
    diag(x)
  }
}
# Convert a linear predictor to the mean scale.
GetMu <- function(Linear_Predictor, Family) {
  family <- .bt_match_family(Family)
  if (is.list(Linear_Predictor)) {
    return(lapply(Linear_Predictor, GetMu, Family = family))
  }
  if (family == "continuous") {
    Linear_Predictor
  } else if (family %in% c("binary", "ordinal")) {
    plogis(Linear_Predictor)
  } else {
    exp(Linear_Predictor)
  }
}
# Convert a linear predictor to the mean scale used in lambda estimation.
GetMu_Lambda <- function(Linear_Predictor, Family) {
  family <- .bt_match_family(Family)
  if (is.list(Linear_Predictor)) {
    return(lapply(Linear_Predictor, GetMu_Lambda, Family = family))
  }
  if (family == "continuous") {
    rep(1, length(Linear_Predictor))
  } else if (family %in% c("binary", "ordinal")) {
    plogis(Linear_Predictor)
  } else {
    exp(Linear_Predictor)
  }
}
# Apply the H transform.
Transform_H <- function(Mu, Family) {
  family <- .bt_match_family(Family)
  if (is.list(Mu)) {
    return(lapply(Mu, Transform_H, Family = family))
  }
  if (family == "continuous") {
    DiagMat(rep(1, length(Mu)))
  } else if (family %in% c("binary", "ordinal")) {
    DiagMat(Mu * (1 - Mu))
  } else {
    DiagMat(Mu)
  }
}
# Approximate matching index.
AppoxMatch <- function(x, y) {
  vapply(seq_along(x), function(i) which.min(abs(x[i] - y)), integer(1))
}
# Remove covariates with all elements in a row or column missing.
RemoveMiss.Fun <- function(X) {
  n <- nrow(X)
  keep.row <- vapply(seq_len(n), function(i) {
    !all(is.na(X[i, , drop = FALSE]))
  }, logical(1))
  row.remove <- which(!keep.row)
  if (length(row.remove) > 0L) {
    X <- X[-row.remove, , drop = FALSE]
  }
  p <- ncol(X)
  keep.col <- vapply(seq_len(p), function(i) {
    !all(is.na(X[, i, drop = FALSE]))
  }, logical(1))
  col.remove <- which(!keep.col)
  if (length(col.remove) > 0L) {
    X <- X[, -col.remove, drop = FALSE]
  }
  list(
    X = X,
    id.remove = if (length(row.remove) > 0L) row.remove else NULL
  )
}
blup.solve <- function(transf.data, membership, sigma, Kmax) {
  lapply(seq_len(Kmax), function(k) {
    pt.k <- membership == k
    XX <- Reduce("+", lapply(which(pt.k), function(j) {
      Xnew <- transf.data[[j]]$Xnew
      t(Xnew) %*% Xnew
    }))
    XY <- Reduce("+", lapply(which(pt.k), function(j) {
      Xnew <- transf.data[[j]]$Xnew
      Ynew <- transf.data[[j]]$Ynew
      t(Xnew) %*% Ynew
    }))
    XZ <- Reduce("+", lapply(which(pt.k), function(j) {
      Xnew <- transf.data[[j]]$Xnew
      Znew <- transf.data[[j]]$Znew
      t(Xnew) %*% Znew
    }))
    ZZ <- Reduce("+", lapply(which(pt.k), function(j) {
      Znew <- transf.data[[j]]$Znew
      t(Znew) %*% Znew
    }))
    ZY <- Reduce("+", lapply(which(pt.k), function(j) {
      Znew <- transf.data[[j]]$Znew
      Ynew <- transf.data[[j]]$Ynew
      t(Znew) %*% Ynew
    }))
    Q <- ZZ + diag(sigma, nrow(ZZ))
    V <- XZ %*% solve(Q, diag(1, nrow(ZZ)))
    A <- XX - V %*% t(XZ)
    b <- XY - V %*% ZY
    fix.eff <- tryCatch(qr.solve(A, b), error = function(ex) NULL)
    if (is.null(fix.eff)) {
      fix.eff <- rep(0, ncol(A))
    }
    rnd.eff <- tryCatch(qr.solve(Q, ZY - t(XZ) %*% fix.eff), error = function(ex) NULL)
    if (is.null(rnd.eff)) {
      rnd.eff <- rep(0, ncol(Q))
    }
    list(fix.eff = fix.eff, rnd.eff = rnd.eff)
  })
}
l2Dist <- function(y1, y2) {
  if (length(y1) != length(y2)) {
    stop("`y1` and `y2` must have the same length.")
  }
  sqrt(mean(unlist(lapply(seq_along(y1), function(i) {
    mean((unlist(y1[[i]]) - unlist(y2[[i]]))^2, na.rm = TRUE)
  })), na.rm = TRUE))
}
line.plot <- function(x, y, ...) {
  mapply(lines, x, y = y, col = "gray", lty = 2)
  invisible(NULL)
}
penBS <- function(d, pen.ord = 2) {
  if (d >= (pen.ord + 1L)) {
    D <- diag(d)
    for (k in seq_len(pen.ord)) {
      D <- diff(D)
    }
    t(D) %*% D
  } else {
    diag(0, d)
  }
}
penBSderiv <- function(d, pen.ord = 2) {
  if (d > 0L) {
    if (d >= (pen.ord + 1L)) {
      pen.matx <- penBS(d, pen.ord)
      cbind(0, rbind(0, pen.matx))
    } else {
      warning("not enough degrees of freedom for differencing penalty matrix: setting penalty to zero\n")
      pen.matx <- diag(1, d + 1L)
      pen.matx[1, 1] <- 0
      pen.matx
    }
  } else {
    0
  }
}
point.plot <- function(x, y, ...) {
  mapply(points, x, y = y, pch = 16)
  invisible(NULL)
}
rho.inv <- function(ni, rho, tol = 1e-2) {
  m <- ni - 1L
  if (m == 0L) {
    0
  } else if (rho < 0 && abs(rho + 1 / m) <= tol) {
    (-1 / m + tol) / (m * tol)
  } else {
    rho / (1 + m * rho)
  }
}
rho.inv.sqrt <- function(ni, rho, tol = 1e-2) {
  m <- ni - 1L
  if (m == 0L) {
    0
  } else {
    if (rho < 0 && abs(rho + 1 / m) <= tol) {
      rho <- -1 / m + tol
    }
    ri <- rho / (1 + m * rho)
    as.numeric(Re(polyroot(c(ri, -2, ni))))[1]
  }
}
sigma.robust <- function(lambda, rho) {
  lambda
}
