\name{simLong}
\alias{simLong}
\title{Simulate longitudinal data}
\description{
  Simulates longitudinal data from models with increasing complexity of
  covariate-time interactions (Ishwaran et al., 2014).
}
\usage{simLong(n,
 ntest = 0,
 N = 5,
 rho = 0.8,
 type = c("corCompSym", "corAR1", "corSymm", "iid"),
 model = c(1, 2, 3),
 phi = 1,
 q = 0,
 ...)}
\arguments{
  \item{n}{Requested training sample size.}
  \item{ntest}{Requested test sample size.}
  \item{N}{Parameter controlling number of time points per subject.}
  \item{rho}{Correlation parameter.}
  \item{type}{Type of correlation matrix.}
  \item{model}{Requested simulation model.}
  \item{phi}{Variance of measurement error.}
  \item{q}{Number of zero-signal variables (i.e., variables unrelated to y).}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Simulates longitudinal data with 3 main effects and a covariate-time
  interaction.  Complexity of the covariate-time interaction
  is specified using the option \code{model}:
  \enumerate{
    \item \emph{\code{model=1}:} Linear covariate-time interaction.
    \item \emph{\code{model=2}:} Quadratic time-quadratic covariate
    interaction.
    \item \emph{\code{model=3}:} Quadratic time-quadratic two-way
    covariate interaction.
  }
  For details see Ishwaran et al. (2014).
}
\value{
  An invisible list with the following components:
  \item{dtaL}{List containing the simulated data in the following order:
    \code{features}, \code{time}, \code{id} and \code{y}.}
  \item{dta}{Simulated data given as a data frame.}
  \item{trn}{Index of \code{id} values identifying the training data.}
  \item{f.true}{Formula of the simulation model.}
}
\author{
  Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran et al. (2014).  Boosting multivariate trees for
  longitudinal data.
}
\keyword{simulation}
