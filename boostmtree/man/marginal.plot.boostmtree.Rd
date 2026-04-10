\name{marginal.plot}
\alias{marginal.plot}
\alias{marginal.plot.boostmtree}
\alias{plot.marginal.plot.boostmtree}
\title{Marginal fitted-response summaries for boostmtree models}
\description{
Create marginal fitted-response summaries for one or more covariates from a
fitted \code{boostmtree} object. The function can draw the resulting curves,
write them to a PDF file, or return the raw and smoothed curve data.
}
\usage{
marginal.plot(object, ...)

\method{marginal.plot}{boostmtree}(
  object,
  M = NULL,
  x.var.names = NULL,
  time.points = NULL,
  subset = NULL,
  prob.class = FALSE,
  response.labels = NULL,
  output = c("plot", "data", "pdf"),
  file = NULL,
  verbose = TRUE,
  ...
)

\method{plot}{marginal.plot.boostmtree}(
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
  \item{x.var.names}{Names of the covariates for which marginal summaries are
  requested. By default all fitted covariates are used.}
  \item{time.points}{Time points at which the fitted curves should be displayed.
  By default the deciles of the observed time values are used. Requested values
  are matched to the nearest observed training times.}
  \item{subset}{Optional vector selecting the fitted subjects to be used in the summary.}
  \item{prob.class}{Logical flag used only for ordinal fits. When \code{TRUE},
  the function works with class probabilities rather than cumulative probabilities.}
  \item{response.labels}{Optional character vector selecting which response
  components to include in the summary. For continuous fits the only available
  label is \code{"response"}. For binary, nominal, and ordinal fits the available
  labels are taken from the fitted object.}
  \item{output}{Either \code{"plot"} to draw the plots on the active graphics
  device, \code{"data"} to return the computed summary object without drawing it,
  or \code{"pdf"} to write the plots to a PDF file.}
  \item{file}{Optional file name used when \code{output = "pdf"} or
  \code{plot(..., output = "pdf")} is requested.}
  \item{verbose}{Should the file location be reported when PDF output is requested?}
  \item{x}{An object returned by \code{marginal.plot}.}
  \item{...}{Further arguments passed to \code{\link{predict.boostmtree}}.}
}
\details{
Marginal plots display the relationship between a covariate and the fitted
response without averaging over modified covariate values. The fitted mean
response is first obtained on a common time grid from
\code{\link{predict.boostmtree}}, and then the association between the chosen
covariate and that fitted response is summarized separately at each requested
time point.

Compared with partial dependence, this is a less adjusted summary because the
other covariates are not perturbed. It is often quicker to compute and is useful
for seeing the broad shape of the fitted longitudinal response surface. When
multiple response components are available, \code{response.labels} can be used
to restrict the display to a chosen subset.
}
\value{
If \code{output = "data"}, an object of class \code{"marginal.plot.boostmtree"}
with components:
\describe{
  \item{\code{data}}{Raw \eqn{x}-versus-fitted-response pairs for each requested
  time point.}
  \item{\code{smooth}}{Smoothed marginal curves used for plotting.}
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
## -------------------------------------------------------------
## Continuous longitudinal marginal summaries.
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
marginal.plot(fit, x.var.names = c("x1", "x2"))

## Return the smoothed marginal curves.
mp <- marginal.plot(
  fit,
  x.var.names = "x2",
  output = "data"
)

str(mp$smooth)

## -------------------------------------------------------------
## Ordinal illustration using class probabilities.
## -------------------------------------------------------------
\donttest{
fit.ord <- boostmtree(
  x = dta$features,
  tm = dta$time,
  id = dta$id,
  y = cut(dta$y, breaks = 3, labels = c("low", "mid", "high")),
  family = "ordinal",
  M = 25,
  verbose = FALSE
)

marginal.plot(
  fit.ord,
  x.var.names = "x1",
  prob.class = TRUE,
  response.labels = c("low", "high")
)
}
}
\keyword{plot}
