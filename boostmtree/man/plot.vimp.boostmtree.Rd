\name{plot.vimp.boostmtree}
\alias{plot.vimp.boostmtree}
\title{Plot Method for Variable Importance Objects}
\description{
Plot variable-importance objects returned by \code{vimp.boostmtree()}.
}

\usage{
\method{plot}{vimp.boostmtree}(
  x,
  show.interaction = TRUE,
  show.time.effect = TRUE,
  output = c("plot", "data", "pdf"),
  file = NULL,
  main = "Variable importance (\%)",
  col = grey(0.80),
  cex.names = 0.8,
  eps = 0.1,
  ...
)
}

\arguments{

  \item{x}{An object of class \code{vimp.boostmtree}.}

  \item{show.interaction}{Logical value indicating whether
    time-interaction effects should be displayed when available.}
  

  \item{show.time.effect}{Logical value indicating whether the time-only
    effect should be reported in the plot margin when available.}
  
  \item{output}{One of \code{"plot"}, \code{"data"}, or
    \code{"pdf"}. The default draws on the active graphics device.}
  

  \item{file}{File name used when \code{output = "pdf"}.}

  \item{main}{Main title for the plot.}

  \item{col}{Bar color.}

  \item{cex.names}{Text magnification for variable names on the x-axis.}

  \item{eps}{Small amount added above the plotted range.}

  \item{...}{Additional arguments passed to \code{barplot()}.}

}

\details{

The default behavior is to draw the variable-importance plot on the
active graphics device.

For longitudinal fits, the plot shows main covariate effects above the x-axis
and time-interaction effects below the x-axis. When available, the time-only
effect is reported in the plot margin.

When \code{output = "data"}, the function returns the plotting data instead of
drawing anything. This can be useful if you want to customize the display.
}
\value{
If \code{output = "plot"} or \code{output = "pdf"}, the function returns the
input object invisibly.

If \code{output = "data"}, the function returns a list containing the plotting
data for each response component.
}
\examples{
set.seed(19)
sim.obj <- simLong(n = 15, n.time = 4, model = 1, family = "continuous")
dta <- sim.obj$data.list

fit <- boostmtree(
  x = dta$features,
  tm = dta$time,
  id = dta$id,
  y = dta$y,
  family = "continuous",
  M = 5,
  cv.flag = TRUE,
  verbose = FALSE,
  control = boostmtree.control(seed = 19)
)

vimp.obj <- vimp.boostmtree(fit, x.names = c("x1", "x2"))
plot(vimp.obj)

plot.data <- plot(vimp.obj, output = "data")
str(plot.data, max.level = 1)
}
