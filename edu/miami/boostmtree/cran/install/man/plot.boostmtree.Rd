\name{plot.boostmtree}
\alias{plot.boostmtree}
\title{Plot Summary Analysis}
\description{
  Plot summary analysis of the boosting analysis.
}
\usage{\method{plot}{boostmtree}(x, ...)}
\arguments{
  \item{x}{An object of class \code{(boostmtree, grow)} or \code{(boostmtree, predict)}.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Plot summary output, including predicted values and residuals.  Also plots
  various parameters against the number of boosting iterations.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran et al. (2014).  Boosting multivariate trees for
  longitudinal data.
}
\keyword{plot}
