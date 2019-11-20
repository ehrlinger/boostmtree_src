simLong <-  function(n = 100,
                     ntest = 0,
                     N = 5,
                     rho = 0.8,
                     type = c("corCompSym", "corAR1", "corSymm", "iid"),
                     model = c(0, 1, 2, 3),
                     family = c("Continuous","Binary"),
                     phi = 1,
                     q = 0,
                     ...)
{
  if(length(family) != 1){
    stop("Specify any one of the family")
  }
  
  if(any(is.na( match(family,c("Continuous","Binary"))))){
    stop("family must be Continuous or Binary")
  }
  type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))
  model <- as.numeric(match.arg(as.character(model), as.character(0:3)))
  dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
    Ni <- round(runif(1, 1, 3 * N))
    type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))
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
    eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
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
    tm <- sample((1:(3 * N))/N, size = Ni, replace = TRUE)
    if (model == 0) {
      l_pred <- 1.5 + 2.5 * x1 - 1.2 * x3 - .6 * x4 + eps
    }
    if (model == 1) {
      l_pred <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * x4 - .65 * tm  * x2   + eps
    }
    if (model == 2) {
      l_pred <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * x4 - .65 * (tm ^ 2) * (x2 ^ 2)   + eps
    }
    if (model == 3) {
      l_pred <- 1.5 + 2.5 * x1 - 1.2 * x3 - .2 * exp(x4) - .65 * (tm ^ 2) * (x2 ^2) * x3  + eps
    }
    mu <- GetMu(Linear_Predictor = l_pred,Family = family)
    if(family == "Continuous"){
      y <- mu
    }
    else {
      y <- rbinom(n = Ni,size = 1,prob = mu)
    }
    cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
          tm, rep(i, Ni), y)
  })))
  d <- q + 4
  colnames(dta) <- c(paste("x", 1:d, sep = ""), "time", "id", "y")
  dtaL <- list(features = dta[, 1:d], time = dta$time, id = dta$id, y = dta$y) 
  if (model == 0) {
    f.true <- "y ~ g( x1 + x3 + x4 )"
  }
  if (model == 1) {
    f.true <- "y ~ g( x1 + x3 + x4 + I(time * x2) )"
  }
  if (model == 2) {
    f.true <- "y ~ g( x1 + x3 + x4 + I(time^2 * x2^2) )"
  }
  if (model == 3) {
    f.true <- "y ~ g( x1 + x3 + exp(x4) + I(time^2 * x2^2 * x3) )"
  }
  trn <- c(1:sum(dta$id <= n))
  return(invisible(list(dtaL = dtaL, dta = dta, trn = trn, f.true = f.true)))
}
