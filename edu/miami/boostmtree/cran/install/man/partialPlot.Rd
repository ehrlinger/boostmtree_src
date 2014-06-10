\name{partialPlot}
\alias{partialPlot}
\title{Partial Plot analysis}
\description{
  Blah blah blah blah blabity blah ...
}
\usage{partialPlot(obj,
 xvar.names,
 tm,
 npts = 25,
 subset,
 plot.it = TRUE,
 ...)}
\arguments{
  \item{obj}{...}
  \item{xvar.names}{...}
  \item{npts}{...}
    \item{tm}{...}
  \item{subset}{...}
  \item{plot.it}{...}
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
