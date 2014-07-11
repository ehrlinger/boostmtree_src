\name{partialPlot}
\alias{partialPlot}
\title{Partial Plot analysis}
\description{
  Partial dependence plot of x against the adjusted y.
}
\usage{partialPlot(obj,
 xvar.names,
 tm,
 npts = 25,
 subset,
 plot.it = TRUE,
 ...)}
\arguments{
  \item{obj}{A boosting object of class \code{(boostmtree, grow)}.}
  \item{xvar.names}{Names of the x-variables to be used.  By default,
    all variables are plotted.}
  \item{tm}{Time points used for the plots of x against y.  By default,
    the deciles of the observed time values are used.}
  \item{npts}{Maximum number of points used for x.  Reduce this value if
    plots are slow.}
  \item{subset}{Vector indicating which rows of the x-data to be used for the
    analysis.  The default is to use the entire data.}
  \item{plot.it}{Should plots be displayed?}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Partial dependence plot (Friedman, 2001) of x values specified by
  \code{xvar.names} against the adjusted y-values over a set of time
  points specified by \code{tm}.  Analysis can be restricted to a subset
  of the data using \code{subset}.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Friedman J.H. Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232, 2001.
}
\examples{
\donttest{
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
colnames(dta) <- c("x1", "x2", "id", "y")

#basic boosting call
boost.grow <- boostmtree(dta[, 1:2], dta[, 3], dta[, 4], dta[, 5], M = 300)

#plot results
partialPlot(boost.grow, "x1")
partialPlot(boost.grow, "x2")

##----------------------------------------------------------------------------
## spirometry data
##----------------------------------------------------------------------------
data(spirometry, package = "_PROJECT_PACKAGE_NAME_")

## boosting call: cubic B-splines with 5 knots 
spr.obj <- boostmtree(spirometry$features, spirometry$time, spirometry$id, spirometry$y,
            M = 200, nu = .05, nknots = 5, d = 3, lambda = -1)

## partial plot of double-lung group at 5 years
dltx <- partialPlot(spr.obj, "AGE", tm = 5, subset=spr.obj$x$DOUBLE==1)

## partial plot of single-lung group at 5 years
sltx <- partialPlot(spr.obj, "AGE", tm = 5, subset=spr.obj$x$DOUBLE==0)

## combine the two plots: we use lowess smoothed values
dltx <- dltx$l.obj[[1]]
sltx <- sltx$l.obj[[1]]
plot(range(c(dltx[, 1], sltx[, 1])), range(c(dltx[, -1], sltx[, -1])),
     xlab = "age", ylab = "predicted y (adjusted)", type = "n")
lines(dltx[, 1], dltx[, -1], lty = 1, lwd = 2, col = "red")
lines(sltx[, 1], sltx[, -1], lty = 1, lwd = 2, col = "blue")
legend("topright", legend = c("DLTx", "SLTx"), lty = 1, fill = c(2,4))
}
}
\keyword{plot}
