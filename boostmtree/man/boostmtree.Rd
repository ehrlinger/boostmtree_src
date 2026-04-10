\name{boostmtree}
\alias{boostmtree}
\title{Boosted multivariate trees for longitudinal data}

\description{

Fit multivariate boosted tree models for repeated-measures outcomes. The
method combines gradient boosting, subject-level tree partitioning, and
smooth time effects to estimate subject-specific mean trajectories or
class probabilities over time.  Supports continuous, binary, nominal,
and ordinal longitudinal data, as well as univariate responses.

}

\usage{
boostmtree(
  x,
  tm = NULL,
  id = NULL,
  y,
  family = c("continuous", "binary", "nominal", "ordinal"),
  y.reference = NULL,
  M = 200,
  nu = 0.05,
  na.action = c("na.omit", "na.impute")[2],
  k = 5,
  mtry = NULL,
  n.knots = 10,
  d = 3,
  pen.ord = 3,
  lambda = NULL,
  rho = NULL,
  lambda.max = 1e6,
  lambda.iter = 2,
  svd.tol = 1e-6,
  verbose = TRUE,
  cv.flag = FALSE,
  eps = 1e-5,
  mod.grad = TRUE,
  nr.iter = 3,
  control = boostmtree.control()
)
}

\arguments{
  
  \item{x}{A data frame or matrix of covariates. For longitudinal fitting,
  rows must be repeated so that the covariate row for subject \eqn{i} appears
  once for each observation time recorded for that subject.}

  \item{tm}{A vector of observation times with one entry per row of \code{x}.
  If \code{tm} and \code{id} are both omitted, the function switches to
  univariate mode and treats each row of \code{x} as a single observation.}

  \item{id}{A subject identifier with one entry per row of \code{x}. Must be
  supplied whenever \code{tm} is supplied.}

  \item{y}{Observed responses with one entry per row of \code{x}.}

  \item{family}{Response type. Choose one of \code{"continuous"},
  \code{"binary"}, \code{"nominal"}, or \code{"ordinal"}.}

  \item{y.reference}{Reference level used when \code{family = "nominal"}.
  If \code{NULL}, the first observed response level is used. Ignored for the
  other families.}

  \item{M}{Number of boosting iterations. Larger values allow a more flexible
  fit, while \code{cv.flag = TRUE} can be used to select a stopping iteration
  automatically.}

  \item{nu}{Boosting step size. Supply either a scalar or a length-two vector.
  The first element is used for the intercept part of the time basis and the
  second for the remaining time-basis coefficients. All values must lie in
  \eqn{(0, 1]}. Smaller values generally require larger \code{M}.}

  \item{na.action}{How to handle missing values in the covariates \code{x}.
  Use \code{"na.omit"} to remove subjects with incomplete covariates, or
  \code{"na.impute"} to keep subjects unless all covariates are missing.
  Missing values in \code{tm}, \code{id}, or \code{y} are not allowed.}

  \item{k}{Requested number of terminal nodes in the tree base learner.}

  \item{mtry}{Number of candidate covariates sampled at each split. If
  \code{NULL}, all covariates are considered.}

  \item{n.knots}{Number of internal knots used for the B-spline time basis.}

  \item{d}{Degree of the B-spline basis. Values smaller than 1 remove the
  smooth time component and leave an intercept-only time design.}

  \item{pen.ord}{Differencing order used in the P-spline penalty. Larger
  values encourage smoother time profiles.}

  \item{lambda}{Smoothing parameter for the time-basis penalty. If
  \code{NULL}, the function estimates \code{lambda} adaptively when the time
  basis has sufficient degrees of freedom. A scalar is recycled across
  subproblems; otherwise the supplied vector must have length \code{n.q}.}

  \item{rho}{Working within-subject correlation parameter. If \code{NULL}, it
  is updated during fitting. A scalar is recycled across subproblems;
  otherwise the supplied vector must have length \code{n.q}.}

  \item{lambda.max}{Upper bound used during adaptive estimation of
  \code{lambda}.}

  \item{lambda.iter}{Number of updates used when \code{lambda} is estimated
  adaptively.}

  \item{svd.tol}{Tolerance used in the singular-value decomposition of the
  penalty matrix.}

  \item{verbose}{Logical; should progress output be printed?}

 \item{cv.flag}{Logical; should out-of-bag cross-validation be tracked to
  select an optimal stopping iteration? When \code{cv.flag = TRUE}, the fit
  must have out-of-bag subjects at each boosting iteration. The default control
  setting \code{control = boostmtree.control()} already provides this through
  \code{bootstrap = "by.root"}. If a no-OOB scheme such as
  \code{bootstrap = "none"} is requested, \code{boostmtree()} switches back
  to \code{"by.root"} with a warning.}

  \item{eps}{Tolerance used when selecting the stopping iteration from the
  cross-validation error path.}

  \item{mod.grad}{Logical; controls the pseudo-response used to grow
    each tree base learner.  During tree growing, the working gradient
    is first computed for each subject and then passed to \code{rfsrc()}
    as the multivariate pseudo-response used to determine tree splits.
    If \code{mod.grad = FALSE}, the inverse of the current working
    covariance matrix is used to compute the gradient, otherwise the
    inverse covariance factor is omitted.  Operationally,
    \code{mod.grad} changes how candidate splits are evaluated, but it
    does \emph{not} remove covariance weighting from the later
    terminal-node coefficient update. After the tree structure has been
    chosen, the node-specific coefficient estimates are still obtained
    from score and Hessian updates that use the inverse working
    covariance matrix.}

  \item{nr.iter}{Number of Newton--Raphson updates used inside the terminal
  node coefficient solver for non-continuous families.}

  \item{control}{A control object created by \code{boostmtree.control()}.}

}

\details{
\code{boostmtree()} is designed for longitudinal data in which subject
\eqn{i} is measured at times \eqn{t_{ij}} and is described by a baseline
covariate vector \eqn{x_i}. The goal is to estimate how the conditional mean
response changes over time and across subjects with different covariate
profiles.

At a high level, the method combines three ideas.

First, it uses \emph{gradient boosting}. Starting from an initial fit, the
algorithm repeatedly adds small updates in the direction that improves the fit.
If \eqn{\eta_i(t)} denotes the current linear predictor, then boosting updates
it according to
\deqn{\eta_i^{(m)}(t) = \eta_i^{(m-1)}(t) + \nu b_m(x_i, t),}
where \eqn{b_m} is the base learner fitted at iteration \eqn{m} and
\eqn{\nu} is the step size.

Second, the base learner is a \emph{multivariate tree}. The tree partitions
subjects into terminal nodes based on their covariates. Within each terminal
node, the time profile is represented by coefficients of a smooth basis in
time. In this way the model can capture nonlinear covariate effects,
interactions, and different trajectory shapes for different groups of
subjects.

Third, the time effect is modeled by a \emph{penalized B-spline}. When
\code{d >= 1}, the observed times are expanded into a spline basis and the
terminal-node coefficients for that basis are penalized to avoid unstable or
overly wiggly trajectories. The smoothing parameter \code{lambda} controls the
strength of that penalty. If \code{lambda} is not supplied, the function can
estimate it adaptively.

The fitted mean is linked to the response scale through
\deqn{g\{\mu_i(t)\} = \eta_i(t),}
where the link \eqn{g} depends on the family:
\itemize{
  \item For \code{"continuous"}, the identity link is used and the method
  models the mean trajectory directly.
  \item For \code{"binary"}, the logit link is used so that the fitted values
  represent estimated probabilities for the non-reference class.
  \item For \code{"nominal"}, the response is decomposed into one-vs-reference
  logit submodels, one for each non-reference class. The final class
  probabilities are then reconstructed from those submodels.
  \item For \code{"ordinal"}, the response is modeled through cumulative logit
  submodels. The final class probabilities are recovered from the fitted
  cumulative probabilities, with a monotonicity correction applied so that the
  fitted cumulative curves remain ordered.
}

For longitudinal data, \code{boostmtree} uses a working covariance model to
describe the dependence among repeated measurements from the same subject after
the mean trajectory has been fitted. For subject \eqn{i} with \eqn{n_i}
observations, and for response component \eqn{q}, the working covariance has
compound-symmetry form

\deqn{
V_{iq} = \phi_q \{ (1 - \rho_q) I_{n_i} + \rho_q J_{n_i} \},
}

where \eqn{I_{n_i}} is the \eqn{n_i \times n_i} identity matrix and
\eqn{J_{n_i}} is the \eqn{n_i \times n_i} matrix of ones.

The roles of \code{phi} and \code{rho} are different:

\itemize{
  \item \code{phi} is a scale parameter. It controls the overall size of the
  residual variation. Larger values of \code{phi} mean that, after accounting
  for the fitted mean structure, the repeated measurements still show greater
  variability.

  \item \code{rho} is a within-subject correlation parameter. It controls how
  strongly repeated measurements from the same subject move together after
  accounting for the fitted mean structure. When \code{rho = 0}, the working
  residuals are uncorrelated within subject. Positive values indicate positive
  within-subject association.
}

Thus, \code{phi} answers the question "how much variation remains?", whereas
\code{rho} answers the question "how is that remaining variation shared across
measurements from the same subject?".

These parameters should be distinguished from \code{lambda}. The parameter
\code{lambda} penalizes roughness in the time-basis coefficients, whereas
\code{phi} and \code{rho} describe the residual covariance structure.

For continuous families there is a single pair \eqn{(\phi, \rho)}. For the
binary, nominal, and ordinal families, the method is built from one or more
working submodels, so the fitted object may contain a separate path
\eqn{(\phi_q, \rho_q)} for each response component \eqn{q}.

Out-of-bag tracking via \code{cv.flag = TRUE} records a standardized error path
and selects an optimal stopping iteration \code{m.opt}. This package uses
\emph{out-of-bag} subjects in two places: for the cross-validation path stored
when \code{cv.flag = TRUE}, and for grow-object variable importance returned by
\code{vimp.boostmtree()}. The default control object uses
\code{bootstrap = "by.root"}, which naturally creates out-of-bag subjects.
By contrast, \code{bootstrap = "none"} produces a deterministic in-bag fit
with no out-of-bag data. That can still be useful for ordinary fitting and for
prediction-object analyses based on a separate test set, but it is not suitable
for out-of-bag CV or grow-object variable importance.

Note that when \code{cv.flag = TRUE}, \code{boostmtree()} checks the requested
resampling rule before fitting. If the requested control object would
yield no out-of-bag subjects, the function either adjusts the resampling
rule to the default out-of-bag scheme or stops with a clear message when
a user-supplied sample matrix is incompatible with OOB CV.

When \code{tm} and \code{id} are omitted, the function switches to univariate
mode. In that case there is no longitudinal time basis, and the method acts as
a boosted tree model for a single response per row.
}


\value{
An object of class \code{c("boostmtree", "grow", ...)} with components:
\describe{
  \item{base.learner}{
    Fitted tree base learners, indexed by boosted subproblem and boosting iteration.
  }

  \item{control}{
    The normalized control object used for fitting.
  }

  \item{cv.flag}{
    Logical; whether cross-validation tracking was requested.
  }

  \item{d}{
    Degree of the time basis used in the fitted model.
  }

  \item{err.rate}{
    Standardized cross-validation error summaries. For single-response fits this is
    typically a matrix with columns \code{"l1"} and \code{"l2"}; for multi-subproblem
    fits it is stored by subproblem. \code{NULL} when \code{cv.flag = FALSE}.
  }

  \item{family}{
    The fitted response family.
  }

  \item{gamma}{
    Terminal-node coefficient summaries by subproblem and boosting iteration.
  }

  \item{gamma.i.list}{
    Cross-validated terminal-node coefficient summaries when
    \code{cv.flag = TRUE}; otherwise \code{NULL}.
  }

  \item{id}{
    The sorted long-format subject identifier.
  }

  \item{id.unique}{
    Unique subject identifiers in subject order.
  }

  \item{k}{
    Requested number of terminal nodes.
  }

  \item{lambda}{
    Path of smoothing-parameter values across boosting iterations and
    subproblems, or \code{NULL} in univariate mode.
  }

  \item{M}{
    Number of boosting iterations.
  }

  \item{m.opt}{
    Selected stopping iteration for each subproblem when
    \code{cv.flag = TRUE}; otherwise \code{NULL}.
  }

  \item{membership}{
    Terminal-node memberships by subproblem and boosting iteration.
  }

  \item{mu}{
    Fitted mean trajectories. For continuous and binary families this is a
    subject-level list of fitted trajectories. For nominal and ordinal families
    this is indexed first by boosted subproblem and then by subject.
  }

  \item{n}{
    Number of subjects.
  }

  \item{n.q}{
    Number of boosted subproblems.
  }

  \item{na.action}{
    Missing-data rule applied to \code{x}.
  }

  \item{ni}{
    Number of observations per subject.
  }

  \item{ntree}{
    Number of trees requested through \code{control}. The current implementation
    requires \code{ntree = 1}.
  }

  \item{nu}{
    Boosting step-size vector used internally.
  }

  \item{oob.available}{
    Logical flag indicating whether the fitted resampling path produced
    out-of-bag subjects at each boosting iteration.
  }

  \item{oob.subject.count}{
    Number of out-of-bag subjects at each boosting iteration. For single-response
    fits this is typically a vector of length \code{M}; for multi-subproblem fits
    it may be stored by subproblem.
  }

  \item{pen.ord}{
    Differencing order used in the P-spline penalty.
  }

  \item{phi}{
    Estimated path of the residual scale parameter. In the working covariance
    model, \code{phi} determines the marginal variance of repeated measurements
    after accounting for the fitted mean structure.
  }

  \item{prob.class}{
    Class probabilities by subject for non-continuous families;
    \code{NULL} for the continuous family.
  }

  \item{q.set}{
    Threshold levels (ordinal) or non-reference levels (binary or nominal)
    used to define the boosted subproblems.
  }

  \item{q.total}{
    Total number of observed response levels for non-continuous families.
  }

  \item{rho}{
    Estimated path of the within-subject correlation parameter. In the working
    covariance model, \code{rho} determines how strongly repeated measurements
    from the same subject are associated after accounting for the fitted mean
    structure. Under compound symmetry, the off-diagonal covariance is
    \eqn{\rho \phi}.
  }

  \item{rmse}{
    Standardized root mean squared error evaluated at \code{m.opt},
    or \code{NULL} when \code{cv.flag = FALSE}.
  }

  \item{time}{
    A list of time vectors, one per subject.
  }

  \item{time.design}{
    Subject-specific time-design matrices.
  }

  \item{time.unique}{
    Sorted unique observation times used to build the time basis.
  }

  \item{univariate}{
    Logical; whether the model was fit in univariate mode.
  }

  \item{x}{
    Subject-level covariate data with one row per subject.
  }

  \item{x.tm}{
    Time-basis design matrix defined on \code{time.unique}.
  }

  \item{x.var.names}{
    Covariate names.
  }

  \item{y}{
    Observed responses split by subject.
  }

  \item{y.levels}{
    Observed response levels for non-continuous families.
    \code{NA} for the continuous family.
  }

  \item{y.mean}{
    Overall response mean used for standardization.
  }

  \item{y.org}{
    Observed responses represented at the boosted-subproblem level.
    For non-continuous families this is the encoded binary or cumulative
    response used by the corresponding submodel.
  }

  \item{y.reference}{
    Reference response level for the nominal family.
    \code{NULL} otherwise.
  }

  \item{y.sd}{
    Overall response standard deviation used for standardization.
  }
}
}



\author{
Amol Pande, Udaya B. Kogalur and Hemant Ishwaran
}

\references{
Friedman J.H. (2001). Greedy function approximation: a gradient boosting
machine. \emph{The Annals of Statistics}, 29(5):1189--1232.

Friedman J.H. (2002). Stochastic gradient boosting.
\emph{Computational Statistics & Data Analysis}, 38(4):367--378.

Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
Blackstone E.H., Ishwaran H. (2017). Boosted multivariate trees for
longitudinal data. \emph{Machine Learning}, 106(2):277--305.

Pande A., Ishwaran H., Blackstone E.H., Rajeswaran J., and Gillanov M.
(2022). Application of gradient boosting in evaluating surgical ablation
for atrial fibrillation. \emph{SN Computer Science}, 3:466.

Pande A., Ishwaran H., and Blackstone E.H. (2022). Boosting for
multivariate longitudinal responses. \emph{SN Computer Science}, 3:186.

Pande A. (2017). \emph{Boosting for longitudinal data}. Ph.D. dissertation,
Miller School of Medicine, University of Miami.
}

\seealso{
  \code{\link{AF}},
  \code{\link{boostmtree.control}},
  \code{\link{marginal.plot.boostmtree}},
  \code{\link{partial.plot.boostmtree}},
  \code{\link{plot.boostmtree}},
  \code{\link{plot.vimp.boostmtree}},
  \code{\link{predict.boostmtree}},
  \code{\link{print.boostmtree}},
  \code{\link{simLong}},
  \code{\link{spirometry}},
  \code{\link{vimp.boostmtree}}
}

\examples{
## -------------------------------------------------------------
## Simulated continuous longitudinal response.
## Use fixed lambda and rho for a fast, example run.
## -------------------------------------------------------------

set.seed(7)
sim.obj <- simLong(n = 20, n.time = 4, rho = 0.5, model = 1,
                   family = "continuous")
dta <- sim.obj$data.list

fit <- boostmtree(
  x = dta$features,
  tm = dta$time,
  id = dta$id,
  y = dta$y,
  family = "continuous",
  M = 10,
  lambda = 0,
  rho = 0,
  verbose = FALSE,
  control = boostmtree.control(seed = 7)
)

print(fit)
plot.data <- plot(fit, output = "data")
str(plot.data[[1]], max.level = 1)

## -------------------------------------------------------------
## Univariate regression: omit tm and id.
## -------------------------------------------------------------

x.uni <- mtcars[, setdiff(names(mtcars), "mpg")]
y.uni <- mtcars$mpg

fit.uni <- boostmtree(
  x = x.uni,
  y = y.uni,
  family = "continuous",
  M = 10,
  lambda = 0,
  verbose = FALSE
)

print(fit.uni)

\donttest{
## -------------------------------------------------------------
## Simulated binary longitudinal response.
## This example uses a fast, deterministic settings.
## -------------------------------------------------------------

set.seed(17)
sim.bin <- simLong(n = 20, n.time = 4, rho = 0.5, model = 2,
                   family = "binary")
dta.bin <- sim.bin$data.list

fit.bin <- boostmtree(
  x = dta.bin$features,
  tm = dta.bin$time,
  id = dta.bin$id,
  y = dta.bin$y,
  family = "binary",
  M = 10,
  lambda = 0,
  rho = 0,
  verbose = FALSE
)

print(fit.bin)
plot(fit.bin)

## -------------------------------------------------------------
## AF data illustration.  Parameters selected for a fast run.
## -------------------------------------------------------------

data(AF, package = "boostmtree")

fit.af <- boostmtree(
  x = AF$features,
  tm = AF$time,
  id = AF$id,
  y = AF$y,
  family = "binary",
  M = 10,
  k = 15,
  nu = .025,
  n.knots = 15,
  verbose = TRUE
)


print(fit.af)


## -------------------------------------------------------------
## Ordinal response illustration.
## A latent continuous response is grouped into three ordered levels.
## -------------------------------------------------------------

set.seed(31)
sim.ord <- simLong(n = 25, n.time = 4, family = "continuous")
dta.ord <- sim.ord$data.list
rank.y <- rank(dta.ord$y, ties.method = "first")
brks <- c(0, floor(length(rank.y) / 3), floor(2 * length(rank.y) / 3), length(rank.y))
y.ord <- cut(
  rank.y,
  breaks = brks,
  labels = c(0, 1, 2),
  include.lowest = TRUE,
  ordered_result = TRUE
)

fit.ord <- boostmtree(
  x = dta.ord$features,
  tm = dta.ord$time,
  id = dta.ord$id,
  y = y.ord,
  family = "ordinal",
  M = 50,
  verbose = FALSE
)

print(fit.ord)

}
}
\keyword{boosting}
