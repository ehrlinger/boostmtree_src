\name{boostmtree}
\alias{boostmtree}
\alias{_PROJECT_PACKAGE_NAME_}
\title{Boosted multivariate trees for longitudinal data}
\description{
  Implements Friedman's gradient descent boosting algorithm for
  longitudinal data using multivariate tree base learners.  A time
  interaction effect is modeled using B-splines and penalized B-splines
  (P-splines) with estimated adaptive smoothing parameter.
}
\usage{
 boostmtree(x,
   tm,
   id,
   y,
   M = 200,
   nu = 0.1,
   K = 3,
   nknots = 10,
   d = 3,
   pen.ord = 2,
   lambda,
   lambda.max = 5e3,
   lambda.iter = 2,
   svd.tol = 1e-6,
   forest.tol = 1e-3,
   verbose = TRUE,
   ...)
}
\arguments{
    \item{x}{Data frame (or matrix) containing the x-values.   Rows
      must be duplicated to match the number of time points for an individual.
      That is, if individual \emph{i} has \emph{n[i]} outcome y-values,
      then there must be \emph{n[i]} duplicate rows of \emph{i}'s x-value.}
    \item{tm}{Vector of time values, one entry for each row in \code{x}.}
    \item{id}{Unique subject identifier, one entry for each row in \code{x}.}
    \item{y}{Observed y-value, one entry for each row in \code{x}.}
    \item{M}{Number of boosting iterations}
    \item{nu}{Boosting regularization parameter.  A value in (0,1].}
    \item{K}{Number of terminal nodes used for the multivariate tree learner.}
    \item{nknots}{Number of knots used for the B-spline for modeling the
    time interaction effect.}
    \item{d}{Degree of the piecewise B-spline polynomial (no time
      effect is fit when d < 1).}
    \item{pen.ord}{Differencing order used to define the penalty with
      increasing values implying greater smoothness.}
    \item{lambda}{Smoothing (penalty) parameter used for B-splines with
      increasing values associated with increasing
      smoothness/penalization. If missing, or non-positive, the value is
      estimated adaptively using a mixed models approach.}
    \item{lambda.max}{Tolerance used for adaptively estimated lambda
      (caps it).  For experts only.}
    \item{lambda.iter}{Number of iterations used to estimate lambda
      (only applies when lambda is not supplied and adaptive smoothing
      is employed).}
    \item{svd.tol}{Tolerance value used in the SVD calculation of the
      penalty matrix.  For experts only.}
    \item{forest.tol}{Tolerance used for forest weighted least squares
      solution.  Experimental and for experts only.}
    \item{verbose}{Should verbose output be printed?}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Each individual has observed y-values over possibly different time
  points, and with possibly differing number of time points.  Given the
  observed y-values, the time points, and the x-value for each
  individual, the mean value of y is estimated using gradient boosting
  in which the gradient is derived from a criterion function involving a
  working variance matrix for y specified as an equicorrelation matrix
  with parameter \emph{rho} multiplied by a variance parameter
  \emph{phi}.  Multivariate trees are used for the base learner and
  weighted least squares are used to solve the line-search.  A
  mixed-model formulation is used to estimate the smoothing parameter
  associated with the B-splines used for the time-interaction effect,
  although the user can manually set the smoothing parameter as well.
}
\value{
     An object of class \code{(boostmtree, grow)} with the following
     components:

     \item{x}{The x-values, but with only one row per individual
       (i.e. duplicated rows are removed). Values sorted on \code{id}.}
     \item{time}{List with each component containing the time
       points for a given individual. Values sorted on \code{id}.}
     \item{id}{Sorted subject identifier.}
     \item{y}{List with each component containing the observed
       y-values for a given individual. Values sorted on \code{id}.}
     \item{ymean}{Overall mean of y-values for all individuals.}
     \item{ysd}{Overall standard deviation of y-values for all individuals.}
     \item{gamma}{List of length \emph{M}, with each component
       containing the boosted tree fitted values.}
     \item{beta}{Matrix of estimated coefficient values for each
       individual.}
     \item{mu}{List with each component containing the estimated mean
       values for an individual.  That is, each component contains the
       estimated time-profile for an individual.}
     \item{lambda}{Smoothing parameter.}
     \item{phi}{Variance parameter.}
     \item{rho}{Correlation parameter.}
     \item{baselearner}{List of length \emph{M} containing the base learners.}
     \item{membership}{List of length \emph{M}, with each component
       containing the terminal node membership for a given boosting
       iteration.}
     \item{D}{Design matrix for time.}
     \item{d}{Degree of the piecewise B-spline polynomial.}
     \item{K}{Number of terminal nodes.}
     \item{M}{Number of boosting iterations.}
     \item{nu}{Boosting regularization parameter.}
     \item{ntree}{Number of trees.}
}
\author{
    Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Friedman J.H. (2001). Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.

  Ishwaran et al. (2014).  Boosting multivariate trees for
  longitudinal data.
}
\seealso{
  \command{\link{partialPlot}}
  \command{\link{plot.boostmtree}},
  \command{\link{predict.boostmtree}},
  \command{\link{print.boostmtree}}
}
\examples{
##------------------------------------------------------------
## Synthetic example
##
##  high correlation, quadratic time with quadratic interaction
##-------------------------------------------------------------
n <- 50
N <- 5
rho <- .80
phi <- 1
dta <- data.frame(do.call("rbind", lapply(1:n, function(i) {
  Ni <- round(runif(1, 1, 3 * N))
  corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
  diag(corr.mat) <- 1
  eps <- sqrt(phi) * t(chol(corr.mat)) \%*\% rnorm(Ni)
  x1 <- rnorm(1)
  x2 <- runif(1, 2, 3)
  x <- c(x1, x2)
  tm <- sample((1:(3 * N))/N, size = Ni, replace = TRUE)
  y <- 1.5 + 2.5 * x1 - .65 * (tm ^ 2) * (x2 ^ 2)  + eps
  cbind(matrix(x, nrow = Ni, ncol = length(x), byrow = TRUE),
        tm, rep(i, Ni), y)
})))
colnames(dta) <- c("x1", "x2", "time", "id", "y")

#basic boosting call (M is set to a small value here purely for illustration)
boost.grow <- boostmtree(dta[, 1:2], dta[, 3], dta[, 4], dta[, 5], M = 100)

#print results
print(boost.grow)

#plot.results
plot(boost.grow)

\donttest{
##----------------------------------------------------------------------------
## spirometry data
##----------------------------------------------------------------------------
data(spirometry, package = "_PROJECT_PACKAGE_NAME_")

## boosting call: B-spline with 3 degree polynomial and 15 knots
spr.obj <- boostmtree(spirometry$features, spirometry$time, spirometry$id, spirometry$y,
                        M = 200, nu = .05, d = 3, nknots = 15)
plot(spr.obj)
}

}
\keyword{boosting}

