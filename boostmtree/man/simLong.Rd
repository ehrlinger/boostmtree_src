\name{simLong}
\alias{simLong}
\title{Simulate longitudinal data for boostmtree examples}
\description{
Generate synthetic longitudinal data for continuous or binary responses.
}
\usage{
simLong(
  n = 100,
  n.test = 0,
  n.time = 5,
  rho = 0.8,
  cor.type = c("cor.comp.sym", "cor.ar1", "cor.symm", "iid"),
  model = c(0, 1, 2, 3),
  family = c("continuous", "binary"),
  phi = 1,
  q = 0,
  ...
)
}
\arguments{
  \item{n}{Number of training subjects.}
  \item{n.test}{Number of additional test subjects to generate.}
  \item{n.time}{Controls the number and scale of repeated observation times.}
  \item{rho}{Working correlation parameter used by the data generator.}
  \item{cor.type}{Correlation structure used for the repeated errors.}
  \item{model}{Signal-generating model. Choose one of \code{0}, \code{1},
  \code{2}, or \code{3}.}
  \item{family}{Response family. Choose \code{"continuous"} or
  \code{"binary"}.}
  \item{phi}{Error-scale parameter.}
  \item{q}{Number of additional noise covariates to append to the baseline
  feature set.}
  \item{...}{Currently ignored.}
}
\value{
A list with components:
\describe{
  \item{data.list}{A list containing \code{features}, \code{time},
  \code{id}, and \code{y}.}
  \item{data}{The full long-format simulated data frame.}
  \item{train.index}{Indices identifying the training portion of
  \code{data}.}
  \item{formula.true}{A character string describing the true signal model.}
}
}
\author{
Hemant Ishwaran, Amol Pande, and Udaya B. Kogalur
}
\seealso{
\code{\link{boostmtree}}
}
\examples{
set.seed(5)
sim.obj <- simLong(n = 10, n.time = 4, family = "continuous")
str(sim.obj$data.list, max.level = 1)
}
