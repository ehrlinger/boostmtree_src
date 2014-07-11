simLong <- function(n = 100,
                     ntest = 0,
                     N = 5,
                     rho = 0.8,
                     type = c("corCompSym", "corAR1", "corSymm", "iid"),
                     model = c(1, 2, 3),
                     phi = 1,
                     q = 0,
                     ...)
{

  ## determine  the type of correlation matrix
  type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))

  ## determine  the requested simulation
  model <- as.numeric(match.arg(as.character(model), as.character(1:3)))

  ## make the data, store as data frame
  dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {

    ## randomly sample the number of time points
    Ni <- round(runif(1, 1, 3 * N))
    type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))

    ## make the correlation matrix (allows different types)
    if (type == "corCompSym") {
      corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
      diag(corr.mat) <- 1
    }
    if (type == "corAR1") {
      corr.mat <- diag(rep(1, Ni))
      if (Ni > 1) {
       for (ii in 1:(Ni - 1)) {
          corr.mat[ii, (ii + 1):Ni] <- rho^(1:(Ni - ii))
        }
        ind <- lower.tri(corr.mat) 
        corr.mat[ind] <- t(corr.mat)[ind]
      }
    }
    if (type == "iid") {
      corr.mat <- diag(rep(1, Ni))
    }

    ## measurement error
    eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
    ## x-variables
    x1 <- rnorm(1)
    x2 <- runif(1, 1, 2)
    x3 <- runif(1, 1, 3)
    x4 <- rnorm(1)
    x <- c(x1, x2, x3, x4)
    p <- length(x)
    if (q > 0) {
      xnoise <- rnorm(q)
      x <- c(x, xnoise)
    }
    
    ## sample the time points
    tm <- sample((1:(3 * N))/N, size = Ni, replace = TRUE)

    ## simulate the y-values
    ## linear time effect
    if (model == 1) {
      y <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * x4 - .65 * tm  * x2   + eps
    }
    ## quadratic-time-xvar
    if (model == 2) {
      y <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * x4 - .65 * (tm ^ 2) * (x2 ^ 2)   + eps
    }
    ## three-way time-xvar
    if (model == 3) {
      y <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * exp(x4) - .65 * (tm ^ 2) * (x2 ^2) * x3  + eps
    }

    ## make the data
    cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
          tm, rep(i, Ni), y)
  })))
  d <- q + 4
  colnames(dta) <- c(paste("x", 1:d, sep = ""), "time", "id", "y")
  
  ## convenient to return the data as a list
  dtaL <- list(features = dta[, 1:d], time = dta$time, id = dta$id, y = dta$y) 

  ## identify the true formula
  if (model == 1) {
    f.true <- "y ~ x1 + x3 + x4 + I(time * x2)"
  }
  if (model == 2) {
    f.true <- "y ~ x1 + x3 + x4 + I(time^2 * x2^2)"
  }
  if (model == 3) {
    f.true <- "y ~ x1 + x3 + exp(x4) + I(time^2 * x2^2 * x3)"
  }
  
  
  ## training data id
  trn <- c(1:sum(dta$id <= n))
  

  ## return the data, training index and formula for the true model
  return(invisible(list(dtaL = dtaL, dta = dta, trn = trn, f.true = f.true)))
  
}
