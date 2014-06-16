##------------------------------------------------------
## internal functions
##------------------------------------------------------
##

## BLUP estimates for random effects (u_k) and fixed effects (alpha_k)
## returned object is a list of lists
blup.solve <- function(transf.data, membership, sigma, Kmax) {

  lapply(1:Kmax, function(k) {
    pt.k <- (membership == k)
    ##sum transformed quantites over the terminal node k
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

    ## use solution on p 19-20 Robinson 1991

    ## .. first we solve for the fixed effects
    Q = ZZ + diag(sigma, nrow(ZZ))
    V = XZ %*% solve(Q, diag(1, nrow(ZZ)))
    A = XX - V %*% t(XZ)
    b = XY - V %*% ZY
    fix.eff <- tryCatch({qr.solve(A, b)}, error = function(ex){NULL})
    if (is.null(fix.eff)) {
      fix.eff <- rep(0, ncol(A))####TBD TBD TBD: what is a good default value?
    }
   
    
    ## .. now solve the random effects
    rnd.eff <- tryCatch({qr.solve(Q, ZY - t(XZ) %*% fix.eff)}, error = function(ex){NULL})
    if (is.null(rnd.eff)) {
      rnd.eff <- rep(0, ncol(Q))####TBD TBD TBD: what is a good default value?
    }

    ## .. test the accuracy of the solution
    ## accr <- sum(abs(XX %*% fix.eff + XZ %*% rnd.eff - XY))

    ## .. return the BLUP solution
    return(list(fix.eff = fix.eff, rnd.eff = rnd.eff))
    
  })
          
}

## constant needed for determining inverse square root
## of equicorrelation matrix
const.sqrt <- function(ni, rho) {
  ri <- rho / (1 + (ni - 1) * rho)
  as.numeric(Re(polyroot(c(ri, -2, ni))))[1]
}


## gls mean-squared error
gls.mse  <- function(f, dta, trn) {
  f <- as.formula(f)
  gls.grow <- tryCatch({gls(f, data = dta[trn, ], correlation = corCompSymm(form = ~ 1 | id))},
                       error = function(ex){NULL})
  if (!is.null(gls.grow)) {
    gls.pred  <- tapply(model.matrix(f, dta[-trn,]) %*% gls.grow$coef,
                        dta[-trn, "id"], function(x) {x})
    y.test <- tapply(dta[-trn, "y"], dta[-trn, "id"], function(x) {x})
    l2Dist(gls.pred, y.test)
  }
  else {
    NA
  }
}


## hidden bootstrap value
is.hidden.bootstrap <-  function (user.option) {

  if (is.null(user.option$bootstrap)) {
    "by.root"
  }
  else {
    as.character(user.option$bootstrap)
  }

}

## hidden ntree value
is.hidden.ntree <-  function (user.option) {

  if (is.null(user.option$ntree)) {
    1
  }
  else {
    max(1, user.option$ntree)
  }

}

## hidden partial plot value
is.hidden.partial <-  function (user.option) {

  if (is.null(user.option$partial)) {
    FALSE
  }
  else {
    TRUE
  }

}

## hidden rho value
is.hidden.rho <-  function (user.option) {

  if (is.null(user.option$rho)) {
    TRUE
  }
  else {
    user.option$rho
  }
  
}


##l1 norm
l1Dist <- function(y1, y2) {
  if (length(y1) != length(y2)) {
    stop("y1 and y2 must have same the length\n")
  }
  mean(unlist(lapply(1:length(y1), function(i) {
    mean(abs(unlist(y1[[i]]) - unlist(y2[[i]])), na.rm = TRUE)
  })), na.rm = TRUE)
}

##l2 norm
l2Dist <- function(y1, y2) {
  if (length(y1) != length(y2)) {
    stop("y1 and y2 must have same the length\n")
  }
  mean(unlist(lapply(1:length(y1), function(i) {
    sqrt(mean((unlist(y1[[i]]) - unlist(y2[[i]]))^2, na.rm = TRUE))
  })), na.rm = TRUE)
}

## modified lowess
lowess.mod <- function(x, y, ...) {
  na.pt <- is.na(x) | is.na(y)
  if (all(na.pt) || sd(y, na.rm = TRUE) == 0) {
    return(list(x = x, y = y))
  }
  else {
    lowess(x[!na.pt], y[!na.pt], ...)
  }
}

## parses boosted tree to determine variable split depth
parse.depth <- function(obj) {
  obj <- stat.split(obj)[[1]]##there is only one tree
  ##determine the depth for each variable
  depth <- unlist(lapply(1:length(obj), function(k) {
    if (!is.null(obj[[k]])) {
      min(obj[[k]][, "dpthID"], na.rm = TRUE)
    }
    else {##did not split
      NA
    }
  }))
  ## assign NA depths (non splitting variables) a maximal depth value
  if (!all(is.na(depth))) {
    ## tree depth
    treeDepth <- max(unlist(lapply(1:length(obj), function(k) {
      if (!is.null(obj[[k]])) {
        max(obj[[k]][, "dpthID"], na.rm = TRUE)
      }
    })), na.rm = TRUE) 
  depth[is.na(depth)] <- treeDepth + 1##maximal depth adds 1 (ONE)
  }
  depth
}

## B-spline penalty
penBS <- function(d, pen.ord = 2) {
  if (d >= (pen.ord + 1)) {
    ## define the differencing matrix
    D <- diag(d)
    for (k in 1:pen.ord) D <- diff(D)
    ## the following yields the penalty matrix P = t(D) %*% D
    t(D) %*% D
  }
  else {
    diag(0, d)
  }
}

## B-spline penalty derivative
penBSderiv <- function(d, pen.ord = 2) {
  if (d >= (pen.ord + 1)) {
    pen.matx <- penBS(d, pen.ord)
    cbind(0, rbind(0, pen.matx))
  }
  else {
    warning("not enough degrees of freedom for differencing penalty matrix: setting penalty to zero\n")
    pen.matx <- diag(1, d + 1)
    pen.matx[1, 1] <- 0
    pen.matx
  }
}


## choose a case at random, plot time profiles for its proximities
plot.profile.prx <- function(obj, col = NULL, rnd.case = NULL, cut = .95, restrictX = TRUE) {
  if (is.null(obj$proximity)) {
    stop("this function requires proximity = TRUE in the predict call")
  }

  ## extract the proximity matix
  prx <- obj$proximity

  ## grow quantities
  time <- obj$boost.obj$time
  time.unq <- sort(unique(unlist(time)))
  DbetaT <- cbind(1, obj$boost.obj$D) %*% t(obj$boost.obj$beta)
  muGrid <- lapply(1:ncol(DbetaT), function(i) {DbetaT[, i]})
  mu <- obj$boost.obj$mu

  ## predict quantities
  time.hat <- obj$time
  muhat <- obj$muhat

  ## proximity calculations
  if (is.null(rnd.case)) {
    rnd.case <- sample(1:nrow(prx), size = 1)
  }
  rnd.prx <- prx[rnd.case, ]
  prx.cut <- quantile(rnd.prx, cut)
  rnd.match <- which(rnd.prx >= prx.cut)
  rnd.prx <- rnd.prx[rnd.match]
  rnd.time <- time.hat[[rnd.case]]
  rnd.mean <- muhat[[rnd.case]]
  
  prx.mu <- c(do.call(cbind, lapply(rnd.match, function(i){muGrid[[i]]}))
                  %*% rnd.prx / sum(rnd.prx))
  prx.which.time <- is.element(time.unq, unlist(lapply(rnd.match, function(i){time[[i]]})))
  prx.time <- time.unq[prx.which.time]
  prx.mu <- prx.mu[prx.which.time] 

  ## custom color
  if (is.null(col)) col <- rep(1, length(rnd.prx))

  ## plot
  plot(supsmu(rnd.time, rnd.mean),
       xlim = if (restrictX) range(rnd.time) else range(c(rnd.time, prx.time)),
       ylim = range(c(rnd.mean, prx.mu, unlist(lapply(rnd.match, function(i){mu[[i]]})))),
       type = "n", xlab = "time", ylab = "mean profile")
  for (i in rnd.match) {
    lines(lowess(time[[i]], mu[[i]]), lty = 2, col = col[i])
    #points(time[[i]], mu[[i]], pch = 16, lty = 2, col = col[i])
  }
  #lines(supsmu(rnd.time, rnd.mean) ,lty = 1, lwd = 2, col = 4)
  points(rnd.time, rnd.mean, pch = 16, cex = 0.25, col = 4)
  lines(lowess(rnd.time, rnd.mean) ,lty = 1, lwd = 2, col = 4)
  points(prx.time, prx.mu, pch = 16, cex = 0.25, col = 1)
  if (restrictX) {
    lines(supsmu(prx.time, prx.mu), lty = 1, lwd = 2, col = 1)
  }
  else {
    lines(lowess(prx.time, prx.mu), lty = 1, lwd = 2, col = 1)
  }
  legend("bottomleft", bty = "n", legend = c(paste("avg prx.", format(mean(rnd.prx, na.rm = TRUE), digits=3))))

  ## return the matched data
  invisible(obj$boost.obj$x[rnd.match,, drop = FALSE])
  
}

