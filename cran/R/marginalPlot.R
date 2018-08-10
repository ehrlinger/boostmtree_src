marginalPlot <- function (obj,
                         xvar.names,
                         tm,
                         subset,
                         plot.it = TRUE,
                         ...)
{
  if (sum(inherits(obj, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)'")
  }
  if (missing(xvar.names)) {
    xvar.names <- colnames(obj$x)
  }
  xvar.names <- intersect(xvar.names, colnames(obj$x))
  if (length(xvar.names) == 0) {
    stop("x-variable names provided do not match original variable names")
  }
  n.xvar <- length(xvar.names)
  tmOrg <- sort(unique(unlist(obj$time)))
  if (missing(tm)) {
    tm.q <- unique(quantile(tmOrg, (1:9)/10, na.rm = TRUE))
    tm.pt <- sapply(tm.q, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  else {
    tm.q <- tm
    tm.pt <- sapply(tm, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  n.tm <- length(tm.pt)
  if (!missing(subset)) {
    obj$x <- obj$x[subset,, drop = FALSE]
  }
  n <- nrow(obj$x)
  if(n.tm == 1){
      muhat <- cbind(matrix(unlist(predict.boostmtree(object = obj,importance = FALSE)$muhat),nrow = n,byrow = TRUE)[,tm.pt])
    }else
      {
      muhat <- matrix(unlist(predict.boostmtree(object = obj,importance = FALSE)$muhat),nrow = n,byrow = TRUE)[,tm.pt]
    }
  lo.obj <- lapply(1:n.xvar, function(nm){
    x <- obj$x[, xvar.names[nm]]
    lo.fit <- lapply(1:n.tm,function(nt){
      fit <- lowess(x,muhat[,nt])
      cbind(fit$x,fit$y)
    })
    names(lo.fit) <- paste("time = ",tm.q,sep="")
    lo.fit
  })
  names(lo.obj) <- xvar.names

  if(plot.it){

    if(n.xvar > 1){
      pdf(file = "MarginalPlot.pdf",width = 10,height = 10)
    }
    for(pp in 1:n.xvar){
      xmin <- min(unlist(lapply(1:n.tm,function(nn){  lo.obj[[pp]][[nn]][,1]   })))
      xmax <- max(unlist(lapply(1:n.tm,function(nn){  lo.obj[[pp]][[nn]][,1]   })))
      ymin <- min(unlist(lapply(1:n.tm,function(nn){  lo.obj[[pp]][[nn]][,2]   })))
      ymax <- max(unlist(lapply(1:n.tm,function(nn){  lo.obj[[pp]][[nn]][,2]   })))

      plot(lo.obj[[pp]][[1]][,1],lo.obj[[pp]][[1]][,2],type = "n",xlim=c(xmin,xmax),ylim=c(ymin,ymax) ,
           xlab = xvar.names[pp],ylab = "Predicted response")
      for(nn in 1:n.tm){
      lines(lo.obj[[pp]][[nn]][,1],lo.obj[[pp]][[nn]][,2],type = "l",col = nn)
      }
    }
    if(n.xvar > 1){
    dev.off()
    print(paste("Plot is stored in the directory:",getwd(),sep=" "))
    }
  }
  return(invisible(lo.obj))
}




