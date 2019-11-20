marginalPlot <- function (object,
                         xvar.names,
                         tm.unq,
                         subset,
                         plot.it = FALSE,
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
  tmOrg <- sort(unique(unlist(object$time)))
  if (missing(tm.unq)) {
    tm.unq <- unique(quantile(tmOrg, (1:9)/10, na.rm = TRUE))
    tm.pt <- sapply(tm.unq, function(tt) {#assign original time values
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
  n <- nrow(object$x)
  if(n.tm == 1){
      muhat <- cbind(matrix(unlist(predict.boostmtree(object = object)$muhat),nrow = n,byrow = TRUE)[,tm.pt])
    }
  else {
      muhat <- matrix(unlist(predict.boostmtree(object = object)$muhat),nrow = n,byrow = TRUE)[,tm.pt]
  }
  
  p.obj <- lapply(1:n.xvar, function(nm){
    x <- object$x[, xvar.names[nm]]
    RawDt <- lapply(1:n.tm,function(nt){
      cbind(x,muhat[,nt])
    })
    names(RawDt) <- paste("time = ",tm.unq,sep="")
    RawDt
  })
  names(p.obj) <- xvar.names  
  
  if(plot.it){
  l.obj <- lapply(1:n.xvar, function(nm){
    x <- object$x[, xvar.names[nm]]
    lo.fit <- lapply(1:n.tm,function(nt){
      fit <- lowess(x,muhat[,nt])
      cbind(fit$x,fit$y)
    })
    names(lo.fit) <- paste("time = ",tm.unq,sep="")
    lo.fit
  })
  names(l.obj) <- xvar.names
  
    if(n.xvar > 1){
      pdf(file = "MarginalPlot.pdf",width = 10,height = 10)
    }
    for(pp in 1:n.xvar){
      xmin <- min(unlist(lapply(1:n.tm,function(nn){  l.obj[[pp]][[nn]][,1]   })))
      xmax <- max(unlist(lapply(1:n.tm,function(nn){  l.obj[[pp]][[nn]][,1]   })))
      ymin <- min(unlist(lapply(1:n.tm,function(nn){  l.obj[[pp]][[nn]][,2]   })))
      ymax <- max(unlist(lapply(1:n.tm,function(nn){  l.obj[[pp]][[nn]][,2]   })))

      plot(l.obj[[pp]][[1]][,1],l.obj[[pp]][[1]][,2],type = "n",xlim=c(xmin,xmax),ylim=c(ymin,ymax) ,
           xlab = xvar.names[pp],ylab = "Predicted response")
      for(nn in 1:n.tm){
      lines(l.obj[[pp]][[nn]][,1],l.obj[[pp]][[nn]][,2],type = "l",col = nn)
      }
    }
    if(n.xvar > 1){
    dev.off()
    print(paste("Plot is stored in the directory:",getwd(),sep=" "))
    }
  }
  return(invisible(list(p.obj = p.obj, l.obj = if(plot.it) l.obj else NULL, time = tmOrg[tm.pt])))
}




