print.boostmtree <- function (x, ...)
{
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("boostmtree", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)' or '(boostmtree, predict)'")
  }
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) == 2) {
    univariate <- length(x$id) == length(unique(x$id))
    cat("model                       :", class(x)[3], "\n")
    cat("fitting mode                :", class(x)[2], "\n")
    cat("Family                      :", x$family, "\n")
    n_levels <- (x$n.Q+1)
    if(x$family == "Nominal" || x$family == "Ordinal"){
    cat("No of levels                :",n_levels, "\n")  
    }
    if (x$ntree > 1) {
      cat("ntree                     :", x$ntree, "\n")
    }
    cat("number of K-terminal nodes  :", x$K, "\n")
    cat("regularization parameter    :", x$nu[1], "\n")
    cat("sample size                 :", nrow(x$x), "\n")
    cat("number of variables         :", ncol(x$x), "\n")
    if (!univariate) {
      cat("number of unique time points:", length(sort(unique(unlist(x$time)))), "\n")
      cat("avg. number of time points  :", round(mean(sapply(x$time, length), na.rm = TRUE), 2), "\n")
      cat("B-spline dimension          :", ncol(x$X.tm), "\n")
      cat("penalization order          :", x$pen.ord, "\n")
    }
    else {
      cat("univariate family           :", TRUE, "\n")
    }
    cat("boosting iterations         :", x$M, "\n")
    if (!is.null(x$err.rate)) {
      if( x$family == "Nominal" || x$family == "Ordinal" ){
       n.Q <- x$n.Q 
      } else
      {
        n.Q <- 1
      }
      optimized_rho <- unlist(lapply(1:n.Q,function(q){
        if(x$family == "Nominal" || x$family == "Ordinal"){
            x$rho[x$Mopt[q],q]
        }else
        {
          x$rho[ x$Mopt[q] ]
        }
      }))
      optimized_phi <- unlist(lapply(1:n.Q,function(q){
        if(x$family == "Nominal" || x$family == "Ordinal"){
           x$phi[x$Mopt[q],q]
        }else
        {
          x$phi[ x$Mopt[q] ]
        }
      }))
      cat("optimized number iterations :", x$Mopt, "\n")      
      if (!univariate) {
        cat("optimized rho               :", round(optimized_rho, 4),  "\n")
        cat("optimized phi               :", round(optimized_phi, 4),  "\n")
      }
        cat("OOB cv RMSE                 :", round(x$rmse, 4), "\n")
    }
  }
  else {
    univariate <- length(x$boost.obj$id) == length(unique(x$boost.obj$id))
    cat("model                       :", class(x)[3], "\n")
    cat("fitting mode                :", class(x)[2], "\n")
    cat("Family                      :", x$family, "\n")
    n_levels <- (x$n.Q+1)
    if(x$family == "Nominal" || x$family == "Ordinal"){
      cat("No of levels                :",n_levels, "\n")  
    }
    cat("sample size                 :", nrow(x$x), "\n")
    cat("number of variables         :", ncol(x$x), "\n")
    if (!univariate) {
      cat("number of unique time points:", length(sort(unique(unlist(x$time)))), "\n")
      cat("avg. number of time points  :", round(mean(sapply(x$time, length), na.rm = TRUE), 2), "\n")
      if (!is.null(x$err.rate)) {
        
        if(x$family == "Nominal" || x$family == "Ordinal"){
          n.Q <- x$n.Q 
        } else
        {
          n.Q <- 1
        }
        optimized_rho <- unlist(lapply(1:n.Q,function(q){
          if(x$family == "Nominal" || x$family == "Ordinal"){
             x$boost.obj$rho[x$Mopt[q],q]
          }else
          {
            x$boost.obj$rho[ x$Mopt[q] ]
          }         
        }))
        optimized_phi <- unlist(lapply(1:n.Q,function(q){
          if(x$family == "Nominal" || x$family == "Ordinal"){
             x$boost.obj$phi[x$Mopt[q],q]
          }else
          {
            x$boost.obj$phi[ x$Mopt[q] ]
          }          
        }))
        cat("optimized number iterations :", x$Mopt, "\n")
        cat("optimized rho               :", round(optimized_rho, 4),  "\n")
        cat("optimized phi               :", round(optimized_phi, 4),  "\n")
        cat("test set RMSE               :", round(x$rmse, 4), "\n")
      }
    }
    else {
      if (!is.null(x$err.rate)) {
        cat("optimized number iterations :", x$Mopt, "\n")
        cat("test set RMSE               :", round(x$rmse, 4), "\n")
      }
      cat("univariate family           :", TRUE, "\n")
    }
  }
}
