boostmtree.diag.matrix <- function(x) {
  if (length(x) == 1L) {
    matrix(x, nrow = 1L, ncol = 1L)
  } else {
    diag(x)
  }
}
boostmtree.l1.dist <- function(y1, y2) {
  if (length(y1) != length(y2)) {
    stop("`y1` and `y2` must have the same length.")
  }
  mean(
    unlist(lapply(seq_along(y1), function(i) {
      mean(abs(unlist(y1[[i]]) - unlist(y2[[i]])), na.rm = TRUE)
    })),
    na.rm = TRUE
  )
}
boostmtree.l2.dist <- function(y1, y2) {
  if (length(y1) != length(y2)) {
    stop("`y1` and `y2` must have the same length.")
  }
  sqrt(
    mean(
      unlist(lapply(seq_along(y1), function(i) {
        mean((unlist(y1[[i]]) - unlist(y2[[i]]))^2, na.rm = TRUE)
      })),
      na.rm = TRUE
    )
  )
}
boostmtree.line.plot <- function(x, y, ...) {
  mapply(lines, x, y = y, col = "gray", lty = 2, SIMPLIFY = FALSE)
  invisible(NULL)
}
boostmtree.point.plot <- function(x, y, ...) {
  mapply(points, x, y = y, pch = 16, SIMPLIFY = FALSE)
  invisible(NULL)
}
boostmtree.lowess <- function(x, y, ...) {
  na.point <- is.na(x) | is.na(y)
  if (all(na.point) || sd(y, na.rm = TRUE) == 0) {
    return(list(x = x, y = y))
  }
  lowess(x[!na.point], y[!na.point], ...)
}
boostmtree.pen.bs <- function(d, pen.ord = 2) {
  if (d >= (pen.ord + 1L)) {
    penalty.matrix <- diag(d)
    for (k in seq_len(pen.ord)) {
      penalty.matrix <- diff(penalty.matrix)
    }
    t(penalty.matrix) %*% penalty.matrix
  } else {
    diag(0, d)
  }
}
boostmtree.pen.bs.deriv <- function(d, pen.ord = 2) {
  if (d > 0) {
    if (d >= (pen.ord + 1L)) {
      penalty.matrix <- boostmtree.pen.bs(d, pen.ord)
      cbind(0, rbind(0, penalty.matrix))
    } else {
      warning(
        "Not enough degrees of freedom for the differencing penalty matrix; ",
        "setting the penalty to zero."
      )
      penalty.matrix <- diag(1, d + 1L)
      penalty.matrix[1, 1] <- 0
      penalty.matrix
    }
  } else {
    0
  }
}
boostmtree.rho.inv <- function(ni, rho, tol = 1e-2) {
  m <- ni - 1L
  if (m == 0L) {
    0
  } else if (rho < 0 && abs(rho + 1 / m) <= tol) {
    (-1 / m + tol) / (m * tol)
  } else {
    rho / (1 + m * rho)
  }
}
boostmtree.rho.inv.sqrt <- function(ni, rho, tol = 1e-2) {
  m <- ni - 1L
  if (m == 0L) {
    0
  } else {
    if (rho < 0 && abs(rho + 1 / m) <= tol) {
      rho <- -1 / m + tol
    }
    rho.inverse <- rho / (1 + m * rho)
    as.numeric(Re(polyroot(c(rho.inverse, -2, ni))))[1]
  }
}
boostmtree.sigma.robust <- function(lambda, rho) {
  lambda
}
boostmtree.blup.solve <- function(transformed.data, membership, sigma, node.count) {
  lapply(seq_len(node.count), function(node.id) {
    node.subject <- membership == node.id
    x.tx <- Reduce("+", lapply(which(node.subject), function(j) {
      x.new <- transformed.data[[j]]$x.new
      t(x.new) %*% x.new
    }))
    x.ty <- Reduce("+", lapply(which(node.subject), function(j) {
      x.new <- transformed.data[[j]]$x.new
      y.new <- transformed.data[[j]]$y.new
      t(x.new) %*% y.new
    }))
    x.tz <- Reduce("+", lapply(which(node.subject), function(j) {
      x.new <- transformed.data[[j]]$x.new
      z.new <- transformed.data[[j]]$z.new
      t(x.new) %*% z.new
    }))
    z.tz <- Reduce("+", lapply(which(node.subject), function(j) {
      z.new <- transformed.data[[j]]$z.new
      t(z.new) %*% z.new
    }))
    z.ty <- Reduce("+", lapply(which(node.subject), function(j) {
      z.new <- transformed.data[[j]]$z.new
      y.new <- transformed.data[[j]]$y.new
      t(z.new) %*% y.new
    }))
    q.matrix <- z.tz + diag(sigma, nrow(z.tz))
    v.matrix <- x.tz %*% solve(q.matrix, diag(1, nrow(q.matrix)))
    a.matrix <- x.tx - v.matrix %*% t(x.tz)
    b.vector <- x.ty - v.matrix %*% z.ty
    fixed.effect <- tryCatch(
      qr.solve(a.matrix, b.vector),
      error = function(ex) NULL
    )
    if (is.null(fixed.effect)) {
      fixed.effect <- rep(0, ncol(a.matrix))
    }
    random.effect <- tryCatch(
      qr.solve(q.matrix, z.ty - t(x.tz) %*% fixed.effect),
      error = function(ex) NULL
    )
    if (is.null(random.effect)) {
      random.effect <- rep(0, ncol(q.matrix))
    }
    list(fixed.effect = fixed.effect, random.effect = random.effect)
  })
}
