print.boostmtree <- function (x, ...)
{

  ## check that object is interpretable
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("boostmtree", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)' or '(boostmtree, predict)'")
  }

  ## grow grint

  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) == 2) {

    cat("model                       :", class(x)[1], class(x)[3], "\n")
    cat("fitting mode                :", class(x)[2], "\n")
    cat("number of K-terminal nodes  :", x$K, "\n")
    cat("regularization parameter    :", x$nu[1], "\n")
    cat("number of unique time points:", length(sort(unique(unlist(x$time)))), "\n")
    cat("avg. number of time points  :", round(mean(sapply(x$time, length), na.rm = TRUE), 2), "\n")
    cat("sample size                 :", nrow(x$x), "\n")
    cat("number of variables         :", ncol(x$x), "\n")
    cat("B-spline dimension          :", ncol(x$D), "\n")
    cat("penalization order          :", ncol(x$pen.ord), "\n")
    
  }

  ## predict print

  else {

    cat("model                       :", class(x)[1], class(x)[3], "\n")
    cat("fitting mode                :", class(x)[2], "\n")
    cat("number of unique time points:", length(sort(unique(unlist(x$time)))), "\n")
    cat("avg. number of time points  :", round(mean(sapply(x$time, length), na.rm = TRUE), 2), "\n")
    cat("sample size                 :", nrow(x$x), "\n")
    cat("number of variables         :", ncol(x$x), "\n")
    if (!is.null(x$err.rate)) {
      cat("optimized number iterations :", x$Mopt, "\n")
      cat("optimized rho               :", round(x$boost.obj$rho[x$Mopt], 4),  "\n")
      cat("optimized phi               :", round(x$boost.obj$phi[x$Mopt], 4),  "\n")
      cat("test set error rate         :", round(x$err.rate[x$Mopt, "l2"], 4), "\n")
    }
    
    
  }
  
  ## reset layout


}
