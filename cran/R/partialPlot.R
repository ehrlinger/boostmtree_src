partialPlot <- function (object,
                         xvar.names,
                         tm.unq,
                         xvar.unq = NULL,
                         npts = 25,
                         subset,
                         conditional.xvars = NULL,
                         conditional.values = NULL,
                         plot.it = FALSE,
                         Variable_Factor = FALSE,
                         ...)
{
  if (sum(inherits(object, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)'")
  }
  if (missing(xvar.names)) {
    xvar.names <- colnames(object$x)
  }
  xvar.names <- intersect(xvar.names, colnames(object$x))
  if (length(xvar.names) == 0) {
    stop("x-variable names provided do not match original variable names")
  }
  n.xvar <- length(xvar.names)
  if ( !is.null(conditional.xvars) && !is.null(conditional.values) ) {
    if (length(conditional.xvars) != length(conditional.values)) {
      stop("conditional x-variable and conditional value vectors are not of same length")
    }
    for (i in 1:length(conditional.xvars)) {
      if( is.factor(object$x[,conditional.xvars[i] ])){
        xuniq <- unique(object$x[,conditional.xvars[i] ])
        if ( !any(xuniq == conditional.values[i] )  ) {
          stop("conditional value for the conditional variable:", conditional.xvars[i], " is not from the original data.")
        }
      }
    }
  }
  tmOrg <- sort(unique(unlist(object$time)))
  if (missing(tm.unq)) {
    tm.q <- unique(quantile(tmOrg, (1:9)/10, na.rm = TRUE))
    tm.pt <- sapply(tm.q, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  else {
    tm.pt <- sapply(tm.unq, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  n.tm <- length(tm.pt)
  if (!missing(subset)) {
    object$x <- object$x[subset,, drop = FALSE]
  }
  if( !is.null(conditional.xvars) && !is.null(conditional.values) ){
    n.cond.xvar <- length(conditional.xvars)
    for(i in 1:n.cond.xvar){
      if(is.factor(object$x[, conditional.xvars[i]  ])){
        object$x[, conditional.xvars[i]  ] <- as.factor(conditional.values[i])
      }else
      {
        object$x[, conditional.xvars[i]  ] <- conditional.values[i]
      }
    }
  }
  p.obj <- lapply(xvar.names, function(nm) {
    x <- object$x[, nm]
    n.x <- length(unique(x))
    if(is.null(xvar.unq)){
      x.unq <- sort(unique(x))[unique(as.integer(seq(1, n.x, length = min(npts, n.x))))]  
    }else
    {
      if(!is.list(xvar.unq)){
        stop("xvar.unq must be a list of length same as xvar.names")
      }
      if( length(xvar.unq) != length(xvar.names) ){
        stop("Length of xvar.unq and xvar.names is different")
      }
      if(!identical(xvar.names,names(xvar.unq))){
        stop("Names of xvar.unq must match with xvar.names")
      }
     x.unq <- xvar.unq[[which(names(xvar.unq) == nm)]] 
    }
    newx <- object$x
    rObj <- t(sapply(x.unq, function(xu) {
      newx[, nm] <- rep(xu, nrow(newx))
      if(Variable_Factor){
        newx[, nm] <- as.factor(newx[, nm])
      }
      mu <- predict(object, x = newx, tm = tmOrg, partial = TRUE, ...)$mu
      mn.x <- colMeans(do.call(rbind, lapply(mu, function(mm) {mm[tm.pt]})))
      c(xu, mn.x)
    }))
    colnames(rObj) <- c("x", paste("y.", 1:length(tm.pt), sep = ""))
    rObj
  })
  names(p.obj) <- xvar.names
  if (plot.it) {
  l.obj <- lapply(p.obj, function(pp) {
    x <- pp[, 1]
    y <- apply(pp[, -1, drop = FALSE], 2, function(yy) {
      lowess(x, yy)$y})
    rObj <- cbind(x, y)
    colnames(rObj) <- c("x", paste("y.", 1:length(tm.pt), sep = ""))
    rObj
  })
  names(l.obj) <- xvar.names
    def.par <- par(no.readonly = TRUE)
    for (k in 1:n.xvar) {
      plot(range(l.obj[[k]][, 1]), range(l.obj[[k]][, -1]), type = "n",
           xlab = xvar.names[k], ylab = "predicted y (adjusted)")
      for (l in 1:n.tm) {
        lines(l.obj[[k]][, 1], l.obj[[k]][, -1, drop = FALSE][, l], type = "l", col = 1)
      }
    }
    par(def.par)
  }
  return(invisible(list(p.obj = p.obj, l.obj = if(plot.it) l.obj else NULL, time = tmOrg[tm.pt])))
}
