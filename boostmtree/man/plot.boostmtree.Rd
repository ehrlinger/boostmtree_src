\name{plot.boostmtree}
\alias{plot.boostmtree}
\title{Plot fitted trajectories and diagnostic summaries for a boostmtree object}
\description{
Display fitted trajectories, error paths, and related diagnostics for objects
produced by \code{boostmtree()} or compatible prediction methods.
}
\usage{
\method{plot}{boostmtree}(
  x,
  use.rmse = TRUE,
  output = c("plot", "data", "pdf"),
  file = NULL,
  width = 10,
  height = 10,
  verbose = TRUE,
  ...
)
}
\arguments{
  \item{x}{An object of class \code{c("boostmtree", "grow", ...)} or
  \code{c("boostmtree", "predict", ...)}.}

  \item{use.rmse}{Logical; if \code{TRUE}, error paths are shown on the
  standardized RMSE scale. If \code{FALSE}, the returned or plotted error path
  is converted to mean squared error.}

  \item{output}{How the result should be returned. The default,
  \code{"plot"}, draws directly on the active graphics device. Use
  \code{"data"} to return the plot data without drawing, or \code{"pdf"} to
  write the plot to a PDF file.}

  \item{file}{File path used when \code{output = "pdf"}. Ignored otherwise.}

  \item{width}{Width of the PDF device in inches when
  \code{output = "pdf"}.}

  \item{height}{Height of the PDF device in inches when
  \code{output = "pdf"}.}

  \item{verbose}{Logical; when \code{TRUE} and \code{output = "pdf"}, print
  the location of the written PDF file.}

  \item{...}{Currently ignored. Included for S3 compatibility.}
}
\details{
The plot method is intended to give a quick visual summary of a fitted
\code{boostmtree} object.

For a grow object, the panels typically show some combination of:
\itemize{
  \item fitted trajectories over time,
  \item residual trajectories,
  \item observed-versus-fitted summaries,
  \item the cross-validation error path when \code{cv.flag = TRUE}, and
  \item the iteration paths for the nuisance parameters \code{rho},
  \code{phi}, and \code{lambda}.
}

For a prediction object, the same plotting framework is used when stored test
errors or variable-importance summaries are available.

Unlike the original implementation, the default behavior is now to draw on the
active graphics device. A PDF file is created only when
\code{output = "pdf"} is requested explicitly.

When \code{output = "data"}, the method returns the numerical information used
to make the plot. This is useful if you want to inspect the fitted curves or
construct your own custom graphics.
}
\value{
If \code{output = "data"}, a list of plot-data objects is returned, one element
per boosted subproblem. Each element can contain components such as
\code{time.by.subject}, \code{fitted.by.subject}, \code{observed.by.subject},
\code{error.path}, \code{rho.path}, \code{phi.path}, \code{lambda.path}, and
\code{m.opt}. Otherwise the same list is returned invisibly after plotting or
writing the PDF file.
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
\code{\link{print.boostmtree}},
\code{\link{simLong}}
}
\examples{
## Fit a small model and draw the default summary plot.
set.seed(23)
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

plot(fit)

## Return the plot data instead of drawing.
plot.data <- plot(fit, output = "data")
names(plot.data[[1]])

\donttest{
## Write the plot to a PDF only when requested explicitly.
out.file <- tempfile(fileext = ".pdf")
plot(fit, output = "pdf", file = out.file, verbose = FALSE)
}
}
\keyword{plot}
