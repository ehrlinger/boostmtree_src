plot.boostmtree <- function (x, ...)
{

  ## check that object is interpretable
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("boostmtree", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)' or '(boostmtree, predict)'")
  }

  ## save the original graphical layout
  def.par <- par(no.readonly = TRUE) 
  
  ## grow plot(s)

  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) == 2) {

    ## extract useful values
    n <- length(x$mu)
    M <- x$M

    layout(rbind(c(1, 4), c(2, 5), c(3, 6)), widths = c(1, 1))

    ## predicted mu versus time
    plot(unlist(x$time), unlist(x$mu), xlab = "time", ylab = "predicted", type = "n")
    for (i in 1:n) {
      lines(x$time[[i]], x$mu[[i]], col = "gray", lty = 2)
    }

    ## residual versus time
    plot(unlist(x$time), unlist(x$y) - unlist(x$mu), xlab = "time", ylab = "residual", type = "n")
    for (i in 1:n) {
      lines(x$time[[i]], x$y[[i]] - x$mu[[i]], col = "gray", lty = 2)
    }

    ## predicted mu versus y
    plot(unlist(x$y), unlist(x$mu), xlab = "y", ylab = "predicted", type = "n")
    for (i in 1:n) {
      lines(x$y[[i]], x$mu[[i]], col = "gray", lty = 2)
    }
    
    ## rho/phi/lambda against M
    plot(1:M, x$rho, ylim = range(lowess.mod(1:M, x$rho)$y),
         xlab = "iterations", ylab = expression(rho), type = "n")
    lines(lowess.mod(1:M, x$rho, f = 5/10))
    plot(1:M, x$phi, ylim = range(lowess.mod(1:M, x$phi)$y),
         xlab = "iterations", ylab = expression(phi), type = "n")
    lines(lowess.mod(1:M, x$phi, f = 5/10))
    plot(1:M, x$lambda, ylim = range(lowess.mod(1:M, x$lambda)$y),
         xlab = "iterations", ylab = expression(lambda), type = "n")
    lines(lowess.mod(1:M, x$lambda, f = 5/10))

    
  }

  ## predict plot(s)

  else {

    ## there are no performance values
    if (is.null(x$err.rate)) {

      ## predicted value versus time
      plot(unlist(x$time), unlist(x$mu), xlab = "time", ylab = "predicted", type = "n")
      for (i in 1:length(x$mu)) {
        lines(x$time[[i]], x$mu[[i]], col = "gray", lty = 2)
      }
      
    }

    
    ## performance plots
    else {
            
      ## extract useful values
      M <- x$boost.obj$M
      Mopt <- x$Mopt

      ## is vimp available?
      if (!is.null(x$vimp)) {
        vimp <- TRUE
        layout(rbind(c(1, 3), c(1, 4), c(2, 5), c(2, 6)), widths = c(1, 1))
      }
      else {
        vimp <- FALSE
        layout(rbind(c(1, 2), c(1, 3), c(1, 4), c(1, 5)), widths = c(1, 1))
      }

        
      ## error rate
      plot(1:M, x$err.rate[, "l2"],
           xlab = "iteration", 
           ylab = "prediction error",
           type = "l", lty = 1)
      abline(v=Mopt, lty = 2, col = 2, lwd = 2)
      
      ## barplot of vimp
      if (vimp) {
        vimp <- 100 * (x$vimp / x$err.rate[Mopt, "l2"])
        barplot(vimp, las = 2, ylab = "vimp (%)", cex.names = 1.0)
      }
      
      ## predicted value versus time
      plot(unlist(x$time), unlist(x$mu), xlab = "time", ylab = "predicted", type = "n")
      for (i in 1:length(x$mu)) {
        lines(x$time[[i]], x$mu[[i]], col = "gray", lty = 2)
      }

      ## rho/phi/lambda against M
      plot(1:M, x$boost.obj$rho, ylim = range(lowess.mod(1:M, x$boost.obj$rho)$y), 
           xlab = "iterations", ylab = expression(rho), type = "n")
      lines(lowess.mod(1:M, x$boost.obj$rho, f = 5/10))
      abline(v=Mopt, lty = 2, col = 2, lwd = 2)
      plot(1:M, x$boost.obj$phi, ylim = range(lowess.mod(1:M, x$boost.obj$phi)$y), 
           xlab = "iterations", ylab = expression(phi), type = "n")
      lines(lowess.mod(1:M, x$boost.obj$phi, f = 5/10))
      abline(v=Mopt, lty = 2, col = 2, lwd = 2)
      plot(1:M, x$boost.obj$lambda, ylim = range(lowess.mod(1:M, x$boost.obj$lambda)$y),
           xlab = "iterations", ylab = expression(lambda), type = "n")
      lines(lowess.mod(1:M, x$boost.obj$lambda, f = 5/10))
      abline(v=Mopt, lty = 2, col = 2, lwd = 2)

      
    }
    
  }
  
  ## reset layout
  par(def.par)

}
