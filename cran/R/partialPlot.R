partialPlot <- function (obj,
                         xvar.names,
                         tm,
                         npts = 25,
                         subset,
                         plot.it = TRUE,
                         ...)
{

  ## check that object is interpretable
  if (sum(inherits(obj, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)'")
  }

  ## determine the desired variables
  if (missing(xvar.names)) {
    xvar.names <- colnames(obj$x)
  }
  xvar.names <- intersect(xvar.names, colnames(obj$x))
  if (length(xvar.names) == 0) {
    stop("x-variable names provided do not match original variable names")
  }
  n.xvar <- length(xvar.names)

  ## what are the desired time values?
  tmOrg <- sort(unique(unlist(obj$time)))
  if (missing(tm)) {
    tm.q <- unique(quantile(tmOrg, (1:9)/10, na.rm = TRUE))
    tm.pt <- sapply(tm.q, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  else {
    tm.pt <- sapply(tm, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  n.tm <- length(tm.pt)

  ## has the user asked to subset the data?
  if (!missing(subset)) {
    obj$x <- obj$x[subset,, drop = FALSE]
  }

  ## iterate over the variables, obtained the partial plot values for the desired time points
  p.obj <- lapply(xvar.names, function(nm) {
    x <- obj$x[, nm]
    n.x <- length(unique(x))
    x.unq <- sort(unique(x))[unique(as.integer(seq(1, n.x, length = min(npts, n.x))))]
    newx <- obj$x
    rObj <- t(sapply(x.unq, function(xu) {
      newx[, nm] <- rep(xu, nrow(newx))
      mu <- predict(obj, x = newx, tm = tmOrg, partial = TRUE)$mu
      mn.x <- colMeans(do.call(rbind, lapply(mu, function(mm) {mm[tm.pt]})))
      c(xu, mn.x)
    }))
    colnames(rObj) <- c("x", paste("y.", 1:length(tm.pt), sep = ""))
    rObj
  })
  names(p.obj) <- xvar.names


  ## create the lowess plot object
  l.obj <- lapply(p.obj, function(pp) {
    x <- pp[, 1]
    y <- apply(pp[, -1, drop = FALSE], 2, function(yy) {
      lowess(x, yy)$y})
    rObj <- cbind(x, y)
    colnames(rObj) <- c("x", paste("y.", 1:length(tm.pt), sep = ""))
    rObj
  })
  names(l.obj) <- xvar.names

  ## plot-it?
  if (plot.it) {
  
    ## save the original graphical layout
    def.par <- par(no.readonly = TRUE) 
    
    ## iterate over the variables, plot the parital plot for the desired time points
    for (k in 1:n.xvar) {
      plot(range(l.obj[[k]][, 1]), range(l.obj[[k]][, -1]), type = "n",
           xlab = xvar.names[k], ylab = "predicted y (adjusted)")
      for (l in 1:n.tm) {
        lines(l.obj[[k]][, 1], l.obj[[k]][, -1, drop = FALSE][, l], type = "l", col = 1)
      }
    }

    ## reset layout
    par(def.par)

  }

  ## return the lowess object for the user to create their own plots
  return(invisible(list(p.obj = p.obj, l.obj = l.obj, time = tmOrg[tm.pt])))
    
}
