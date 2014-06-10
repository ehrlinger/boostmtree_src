\name{plot.boostmtree}
\alias{plot.boostmtree}
\title{Plot Error Rate and Variable Importance from a RF-SRC analysis}
\description{
  Plot out-of-bag (OOB) error rates and variable importance (VIMP)
  from a RF-SRC analysis. This is the default plot method for the package.
}
\usage{\method{plot}{boostmtree}(x, ...)}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc, predict)}.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Plot cumulative OOB error rates as a function of number of trees.
  Plot variable importance (VIMP) if available.
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
}
\seealso{
  \command{\link{predict.rfsrc}},
  \command{\link{rfsrc}}
}
\examples{
}
\keyword{plot}
