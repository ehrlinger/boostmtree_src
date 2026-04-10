\name{partial.plot}
\alias{partial.plot}
\alias{partial.plot.boostmtree}
\alias{plot.partial.plot.boostmtree}
\title{Partial-dependence summaries for fitted boostmtree models}
\description{
Create partial-dependence summaries for one or more covariates from a fitted
\code{boostmtree} object. The function can either draw the plots immediately,
write them to a PDF file, or return the computed curve data for further use.
}
\usage{
partial.plot(object, ...)

\method{partial.plot}{boostmtree}(
  object,
  M = NULL,
  x.var.names = NULL,
  time.points = NULL,
  x.var.values = NULL,
  n.points = 25,
  subset = NULL,
  prob.class = FALSE,
  response.labels = NULL,
  conditional.x.var.names = NULL,
  conditional.values = NULL,
  output = c("plot", "data", "pdf"),
  file = NULL,
  verbose = TRUE,
  use.cv.flag = NULL,
  ...
)

\method{plot}{partial.plot.boostmtree}(
  x,
  output = c("plot", "pdf"),
  file = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
  \item{object}{A fitted object of class \code{(boostmtree, grow)}.}
  \item{M}{Number of boosting iterations to use. By default the function uses
  the fitted stopping value when available and otherwise uses the full fitted path.}
  \item{x.var.names}{Names of the covariates for which partial-dependence
  summaries are requested. By default all fitted covariates are used.}
  \item{time.points}{Time points at which the fitted curves should be displayed.
  By default the deciles of the observed time values are used. Requested values
  are matched to the nearest observed training times.}
  \item{x.var.values}{Optional named list giving the covariate values at which
  each partial-dependence curve should be evaluated. When omitted, the function
  uses evenly spaced observed values.}
  \item{n.points}{Maximum number of covariate values used for each partial-dependence
  curve when \code{x.var.values} is not supplied.}
  \item{subset}{Optional vector selecting the fitted subjects to be used when
  averaging the partial-dependence curves.}
  \item{prob.class}{Logical flag used only for ordinal fits. When \code{TRUE},
  the function works with class probabilities rather than cumulative probabilities.}
  \item{response.labels}{Optional character vector selecting which response
  components to include in the summary. For continuous fits the only available
  label is \code{"response"}. For binary, nominal, and ordinal fits the available
  labels are taken from the fitted object.}
  \item{conditional.x.var.names}{Optional names of additional covariates whose
  values should be held fixed while the partial-dependence curves are computed.}
  \item{conditional.values}{Values corresponding to \code{conditional.x.var.names}.}
  \item{output}{Either \code{"plot"} to draw the plots on the active graphics
  device, \code{"data"} to return the computed curve object without drawing it,
  or \code{"pdf"} to write the plots to a PDF file.}
  \item{file}{Optional file name used when \code{output = "pdf"} or
  \code{plot(..., output = "pdf")} is requested.}
  \item{verbose}{Should the file location be reported when PDF output is requested?}
  \item{use.cv.flag}{Should the partial-dependence curves be based on the
  OOB/CV coefficient path? The default is to use the fitted object setting.
  This option is only meaningful when the model was fitted with
  \code{cv.flag = TRUE}; for partial plots it is supported only when all fitted
  subjects are used in their original order.}
  \item{x}{An object returned by \code{partial.plot}.}
  \item{...}{Further arguments passed to \code{\link{predict.boostmtree}}.}
}
\details{
Partial-dependence plots show how the fitted mean response changes as one
covariate is varied over a grid of values while the remaining covariates are
held fixed at their observed values. In a longitudinal model this produces a
family of curves, one for each requested time point.

For each chosen covariate value, the function constructs a modified subject-level
covariate data set, obtains fitted values from \code{\link{predict.boostmtree}},
and then averages those fitted values over subjects at the selected times.
This gives an adjusted view of the fitted effect of the chosen covariate, in the
sense of Friedman (2001), while retaining the time-varying structure of the
longitudinal fit.

The returned object stores the computed curve data and can be plotted later with
\code{plot()}. For families with multiple boosted subproblems or multiple class
probabilities, \code{response.labels} can be used to restrict attention to a
subset of response components.
}
\value{
If \code{output = "data"}, an object of class \code{"partial.plot.boostmtree"}
with components:
\describe{
  \item{\code{curves}}{Partial-dependence curves. For continuous and binary
  fits this is a named list indexed by covariate name. For nominal and ordinal
  fits it is a list indexed first by response component and then by covariate.}
  \item{\code{time.points}}{The observed time points used in the summary.}
  \item{\code{x.var.names}}{Covariate names included in the summary.}
  \item{\code{response.labels}}{Selected labels for the response component(s)
  included in the summary.}
}
If \code{output = "plot"} or \code{"pdf"}, the same object is returned
invisibly after plotting.
}
\references{
Friedman, J. H. (2001). Greedy function approximation: a gradient boosting
machine. \emph{Annals of Statistics}, 29, 1189--1232.
}
\examples{
\donttest{
## -------------------------------------------------------------
## Continuous longitudinal partial dependence.
## -------------------------------------------------------------
set.seed(19)
sim.obj <- simLong(n = 40, n.time = 4, model = 2, family = "continuous")
dta <- sim.obj$data.list

fit <- boostmtree(
  x = dta$features,
  tm = dta$time,
  id = dta$id,
  y = dta$y,
  family = "continuous",
  M = 50,
  verbose = FALSE
)

## Draw directly on the active device.
partial.plot(fit, x.var.names = c("x1", "x2"))

## Return the underlying curve data.
pp <- partial.plot(
  fit,
  x.var.names = "x2",
  output = "data"
)

str(pp$curves)

## -------------------------------------------------------------
## Ordinal illustration using class probabilities.
## -------------------------------------------------------------

fit.ord <- boostmtree(
  x = dta$features,
  tm = dta$time,
  id = dta$id,
  y = cut(dta$y, breaks = 3, labels = c("low", "mid", "high")),
  family = "ordinal",
  M = 25,
  verbose = FALSE
)

partial.plot(
  fit.ord,
  x.var.names = "x1",
  prob.class = TRUE,
  response.labels = c("low", "high")
)
}
}
\keyword{plot}
