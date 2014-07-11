\name{predict.boostmtree}
\alias{predict.boostmtree}
\title{Prediction for Boosted multivariate trees for longitudinal data.}
\description{
   Obtain predicted values.  Also returns test-set performance
   if the test data contains y-outcomes.
}
\usage{\method{predict}{boostmtree}(object,
 x,
 tm,
 id,
 y,
 M,
 importance = TRUE,
 proximity = FALSE,
 verbose = TRUE,
 eps = 1e-4,
 ...)
}
\arguments{
  \item{object}{A boosting object of class \code{(boostmtree, grow)}.}
  \item{x}{Data frame (or matrix) containing test set x-values.  Rows
    must be duplicated to match the number of time points for an
    individual. If missing, the training x values are used and
    \code{tm}, \code{id} and \code{y} are not required.}
  \item{tm}{Time values for each test set individual with one entry for
    each row of \code{x}.  Optional, but if missing, the set of unique time
    values from the training values are used for each individual and no
    test-set performance values are returned.}
  \item{id}{Unique subject identifier, one entry for each row in
    \code{x}.  Optional, but if missing, each individual is assumed to
    have a full time-profile specified by the unique time values from
    the training data.  No test-set performance values are returned.}
  \item{y}{Test set y-values, with one entry for each row  in \code{x}.}
  \item{M}{Fixed value for the boosting step number.  Leave this empty 
    to determine the optimized value obtained by minimizing test-set error.}
  \item{importance}{Calculate importance values for the x-variables?}
  \item{proximity}{Calculate proximity of the test data to the training data?}
  \item{verbose}{Should verbose output be printed?}
  \item{eps}{Tolerance value used for determing the optimal \code{M}.
    For experts only.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  The predicted time profile and performance values are obtained for
  test data from the boosted object grown on the training data.
}
\value{
     An object of class \code{(boostmtree, predict)}, which is a list with the
     following components:
     \item{obj}{The original boosting object.}
     \item{x}{The test x-values, but with only one row per individual
       (i.e. duplicated rows are removed).}
     \item{time}{List with each component containing the time
       points for a given test individual.}
     \item{time.unq}{Vector of the unique time points.}
     \item{y}{List containing the test y-values.}
     \item{beta}{Estimated beta coefficients at the optimized
       \code{M}.}
     \item{mu}{Estimated mean profile at the optimized \code{M}.}
     \item{err.rate}{Test set l1 and l2 (mean-squared error, mse) values.}
     \item{mse}{Test set mean-squared value at the optimized \code{M}.}
     \item{vimp}{Variable importance values.}
     \item{proximity}{Proximity matrix.}
     \item{Mopt}{The optimized \code{M}.}
}
\author{
    Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran et al. (2014).  Boosting multivariate trees for
  longitudinal data.
}
\seealso{
  \command{\link{plot.boostmtree}},
  \command{\link{print.boostmtree}}
}
\examples{
\donttest{
##------------------------------------------------------------
## Synthetic example
##
##  high correlation, quadratic time with quadratic interaction
##  noisy variables
##
##  calculate test set error under adaptive rho and mispecified rho
##  compare to GLS (generalized least-squares)
##----------------------------------------------------------------------------
M <- 250
n <- 100
ntest <- 100
N <- 5
rho <- .80
phi <- 1
p <- 5
dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
  Ni <- round(runif(1, 1, 3 * N))
  corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
  diag(corr.mat) <- 1
  eps <- sqrt(phi) * t(chol(corr.mat)) \%*\% rnorm(Ni)
  x1 <- rnorm(1)
  x2 <- runif(1, 2, 3)
  xnoise <- rnorm(p - 2)
  x <- c(x1, x2, xnoise)
  tm <- sample((1:(3 * N))/N, size = Ni, replace = TRUE)
  y <- 1.5 + 2.5 * x1 - .65 * (tm ^ 2) * (x2 ^ 2)  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        tm, rep(i, Ni), y)
})))
colnames(dta) <- c(paste("x", 1:p, sep = ""), "time", "id", "y")
trn <- c(1:sum(dta$id <= n))

##adaptive rho
boost.grow.1 <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3],
             M = M, lambda = 1000)
boost.pred.1 <- predict(boost.grow.1,
             dta[-trn, 1:p], dta[-trn, p+1], dta[-trn, p+2], dta[-trn, p+3])

##fixed rho = 0
boost.grow.2 <- boostmtree(dta[trn, 1:p], dta[trn, p+1], dta[trn, p+2], dta[trn, p+3], 
             rho = 0, M = M, lambda = 1000)
boost.pred.2 <- predict(boost.grow.2,
             dta[-trn, 1:p], dta[-trn, p+1], dta[-trn, p+2], dta[-trn, p+3])

##mse values
cat("true model        :", boostmtree:::gls.mse("y~x1+I(x2^2)*I(time^2)",dta,trn),"\n")
cat("linear GLS        :", boostmtree:::gls.mse("y~x1+x2+x1*time+x2*time",dta,trn),"\n")
cat("boostmtree        :", boost.pred.1$mse,"\n")
cat("boostmtree (rho=0):", boost.pred.2$mse,"\n")

##predicted value plot
plot(boost.pred.1)
}
}
\keyword{predict}
\keyword{boosting}

