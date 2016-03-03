plot.boostmtree <- function (x, ...)
{

  ## check that object is interpretable
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("boostmtree", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)' or '(boostmtree, predict)'")
  }


  ###########################################################################
  ##
  ## grow plot(s)
  ##
  ###########################################################################
  
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) == 2) {

    ## save the original graphical layout
    def.par <- par(no.readonly = TRUE) 

    ## extract useful values
    n <- length(x$mu)
    M <- x$M
    univariate <- length(x$id) == length(unique(x$id))

    ## plot layout depends on availability of training vimp and if this is a univariate problem
    if (!univariate) {
      if (is.null(x$vimp)) {
        layout(rbind(c(1, 4), c(2, 5), c(3, 6)), widths = c(1, 1))
      }
      else {
        layout(rbind(c(1, 3), c(2, 4), c(2, 5)), widths = c(1, 1))
      }
    }
    else {
      if (!is.null(x$vimp)) {
        layout(rbind(c(1, 2)), widths = c(1, 1))
      }
    }

    ## predicted mu versus time
    if (!univariate) {
      plot(unlist(x$time), unlist(x$mu), xlab = "time", ylab = "predicted", type = "n")
      for (i in 1:n) {
        lines(x$time[[i]], x$mu[[i]], col = "gray", lty = 2)
      }
    }

    ## residuals or vimp
    if (!univariate) {
      if (is.null(x$vimp)) {
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
      }
      else {#vimp
        vimp <- 100 * x$vimp
        barplot(vimp, las = 2, ylab = "vimp (%)", cex.names = 1.0)
      }
    }
    else {
      if (is.null(x$vimp)) {
        ## predicted mu versus y
        plot(unlist(x$y), unlist(x$mu), xlab = "y", ylab = "predicted", type = "n")
        for (i in 1:n) {
          points(x$y[[i]], x$mu[[i]], pch = 16)
        }
        abline(0, 1, col = "gray", lty = 2)
      }
      else {#vimp
        vimp <- 100 * x$vimp
        barplot(vimp, las = 2, ylab = "vimp (%)", cex.names = 1.0)
      }
    }
      
    ## rho/phi/lambda against M
    if (!univariate) {
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

    ## reset layout
    par(def.par)
    
  }

  ###########################################################################
  ##
  ## predict plot(s)
  ##
  ###########################################################################

  else {

    ## do we have a univariate setting?
    univariate <- length(x$boost.obj$id) == length(unique(x$boost.obj$id))

    ## save the original graphical layout
    if (!(univariate && is.null(x$err.rate))) {
      def.par <- par(no.readonly = TRUE) 
    }
    
    ## there are no performance values
    if (!univariate && is.null(x$err.rate)) {

      ## predicted value versus time
      plot(unlist(x$time), unlist(x$mu), xlab = "time", ylab = "predicted", type = "n")
      for (i in 1:length(x$mu)) {
          lines(x$time[[i]], x$mu[[i]], col = "gray", lty = 2)
      }
      
    }

    
    ## performance plots
    else if (!is.null(x$err.rate)) {

      ## extract useful values
      M <- x$boost.obj$M
      Mopt <- x$Mopt

      ## is vimp available?
      if (!univariate) {
        
        if (!is.null(x$vimp)) {
          layout(rbind(c(1, 3), c(1, 4), c(2, 5), c(2, 6)), widths = c(1, 1))
        }
        else {
          layout(rbind(c(1, 2), c(1, 3), c(1, 4), c(1, 5)), widths = c(1, 1))
        }
      }
      else {
        if (!is.null(x$vimp)) {
          layout(rbind(c(1, 2)), widths = c(1, 1))
        }
      }
        
      ## error rate
      plot(1:M, x$err.rate[, "l2"],
           xlab = "iteration", 
           ylab = "RMSE prediction error",
           type = "l", lty = 1)
      abline(v = Mopt, lty = 2, col = 2, lwd = 2)
      
      ## barplot of vimp
      if (!is.null(x$vimp)) {
        vimp <- 100 * (x$vimp / x$err.rate[Mopt, "l2"])
        barplot(vimp, las = 2, ylab = "vimp (%)", cex.names = 1.0)
      }

      if (!univariate) {
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
    if (!(univariate && is.null(x$err.rate))) {
      par(def.par)
    }
  
    
  }
  

}
