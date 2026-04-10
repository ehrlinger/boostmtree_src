\name{print.boostmtree}
\alias{print.boostmtree}
\title{Print a summary of a boostmtree fit or prediction object}
\description{
Print a compact summary of a fitted \code{boostmtree} model or a compatible
prediction object.
}
\usage{
\method{print}{boostmtree}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{c("boostmtree", "grow", ...)} or
  \code{c("boostmtree", "predict", ...)}.}
  \item{...}{Currently ignored. Included for S3 compatibility.}
}
\details{
The printed summary is intended as a quick check that the model you fit is the
model you meant to fit. It reports the response family, the requested number of
terminal nodes, the sample size, and the structure of the repeated-measures
design.

When cross-validation or test-set error information is available, the summary
also reports the selected stopping iteration and the corresponding RMSE. For
longitudinal fits it can additionally report the estimated working-correlation
and working-variance parameters at the selected stopping iteration.
}
\value{
The input object, returned invisibly.
}
\references{
Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
Blackstone E.H., Ishwaran H. (2017). Boosted multivariate trees for
longitudinal data. \emph{Machine Learning}, 106(2):277--305.

Pande A., Ishwaran H., Blackstone E.H., Rajeswaran J., and Gillanov M.
(2022). Application of gradient boosting in evaluating surgical ablation
for atrial fibrillation. \emph{SN Computer Science}, 3:466.

Pande A., Ishwaran H., and Blackstone E.H. (2022). Boosting for
multivariate longitudinal responses. \emph{SN Computer Science}, 3:186.
}
\seealso{
\code{\link{boostmtree}},
\code{\link{plot.boostmtree}},
\code{\link{simLong}}
}
\examples{
## Print a basic longitudinal fit.
set.seed(11)
sim.obj <- simLong(n = 15, n.time = 4, family = "continuous")
dta <- sim.obj$data.list

fit <- boostmtree(
  x = dta$features,
  tm = dta$time,
  id = dta$id,
  y = dta$y,
  family = "continuous",
  M = 8,
  verbose = FALSE
)

print(fit)

## Print a fit that includes out-of-bag stopping information.
set.seed(12)
sim.bin <- simLong(n = 15, n.time = 4, family = "binary")
dta.bin <- sim.bin$data.list

fit.bin <- boostmtree(
  x = dta.bin$features,
  tm = dta.bin$time,
  id = dta.bin$id,
  y = dta.bin$y,
  family = "binary",
  M = 10,
  cv.flag = TRUE,
  verbose = FALSE
)

print(fit.bin)
}
\keyword{print}
