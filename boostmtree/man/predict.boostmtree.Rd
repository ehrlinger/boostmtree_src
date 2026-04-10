\name{predict.boostmtree}
\alias{predict.boostmtree}
\title{Predict longitudinal trajectories from a fitted boostmtree model}
\description{
Generate fitted or predicted trajectories from a \code{boostmtree} fit.
The method can return predictions for the original training subjects,
for a new longitudinal test set, or for new subjects evaluated on a
common time grid.
}
\usage{
\method{predict}{boostmtree}(
  object,
  x,
  tm,
  id,
  y,
  M = NULL,
  eps = 1e-5,
  use.cv.flag = FALSE,
  partial = FALSE,
  ...
)
}
\arguments{

  \item{object}{A fitted object returned by \code{\link{boostmtree}}.}

  \item{x}{New covariate values. Supply either one row per subject or a
    long-format data frame with one row per observation. If omitted, the
    training subjects stored in \code{object} are used.}

  \item{tm}{Observation times for the prediction data. When \code{x},
    \code{tm}, and \code{id} are supplied together, predictions are made at
    the observed times in the new data. When \code{x} is supplied without
    \code{tm} and \code{id}, each row of \code{x} is treated as a new subject
    and predictions are returned on the training time grid.}

  \item{id}{Subject identifier for long-format prediction data.}

  \item{y}{Observed responses for the prediction data. When supplied with
    long-format prediction data, the method computes standardized test-set
    error summaries and, if \code{M = NULL}, chooses the stopping iteration by
    minimizing the test-set error path.}

  \item{M}{Optional fixed boosting iteration used for prediction. If
    \code{NULL} and \code{y} is supplied, the stopping iteration is selected by
    minimizing the test-set error path. If \code{NULL} and \code{y} is not
    supplied, the stored \code{m.opt} from the fitted object is used when
    available; otherwise the full fitted path is used.}

  \item{eps}{Tolerance used when selecting \code{m.opt} from the test-set
    error path.}
  
  \item{use.cv.flag}{Logical; should the prediction use the stored
    out-of-bag coefficient estimates from a model fit with
    \code{cv.flag = TRUE}? This option is intended for predictions on the
    original training subjects and is ignored for new prediction data.}
  
  \item{partial}{Logical; if \code{TRUE}, interpret \code{x} as one row per
    subject and \code{tm} as a common time grid at which each subject should be
    evaluated.}

\item{...}{Currently ignored. Included for S3 compatibility.}

}
\details{

For longitudinal prediction on a new data set, supply \code{x}, \code{tm},
and \code{id} in long format, with the same subject-level covariates used in
the training fit. The function collapses the covariates to one row per subject
and then uses the stored terminal-node coefficient path to reconstruct the
predicted trajectory.

There are three common use cases.

First, if \code{x} is omitted, the method returns fitted trajectories for the
training subjects. This is useful for inspecting the stored fit or for
producing training-set plots.

Second, if \code{x} is supplied without \code{tm} and \code{id}, each row of
\code{x} is treated as a new subject and predictions are returned on the
training time grid stored in \code{object$time.unique}. This is convenient
when the user wants an entire fitted profile for new subjects.

Third, if \code{partial = TRUE}, \code{x} must contain one row per subject and
\code{tm} supplies a common prediction grid. This returns fitted profiles on
that user-specified grid. The method uses the fitted time basis from the
training model, so predictions are meaningful when the supplied grid remains
within or close to the observed training-time range.

For binary and nominal families, prediction proceeds through one-vs-reference
submodels; for ordinal families, prediction proceeds through cumulative
submodels followed by a monotonicity correction across thresholds. The returned
\code{mu} component stores the predicted mean path for each boosted subproblem,
and \code{prob.class} converts these to class probabilities for the original
response scale.

When \code{y} is supplied for the prediction data, the method computes
standardized RMSE along the prediction path.  If \code{M = NULL}, the
method selects \code{m.opt} by minimizing the test-set RMSE, using the
same tolerance rule as in model fitting. }

\value{
An object of class \code{c("boostmtree", "predict", ...)} with components:
\describe{
  \item{base.learner}{Stored tree learners from the fitted object.}

  \item{boost.obj}{The fitted training object with large internal fitting
  components removed.}

  \item{df.time.design}{Number of columns in the time-design matrices.}

  \item{err.rate}{Standardized test-set error summaries. For single-response
  fits this is typically a matrix with columns \code{"l1"} and \code{"l2"};
  for multi-subproblem fits it is stored by subproblem. \code{NULL} when
  prediction responses \code{y} are not supplied.}

  \item{family}{The fitted response family.}

  \item{gamma}{Stored terminal-node coefficient summaries used for prediction.}

  \item{id}{Long-format subject identifier corresponding to the supplied
  prediction data.}

  \item{id.unique}{Unique subject identifiers in subject order.}

  \item{k}{Number of terminal nodes requested during fitting.}

  \item{m.opt}{Selected stopping iteration for each boosted subproblem.}

  \item{membership}{Predicted terminal-node memberships for each subproblem
  and boosting iteration.}

  \item{mu}{Predicted mean trajectories at the time points requested.
     If \code{x}, \code{tm}, and \code{id} are supplied, then \code{mu}
     is evaluated at the supplied subject-specific times. If only
     \code{x} is supplied, then each new subject is predicted on the
     fitted training time grid, so \code{mu} is already a full profile
     on that grid. If \code{partial = TRUE}, then \code{mu} is evaluated
     on the user-supplied common grid \code{tm}.  For continuous and
     binary families this is a subject-level list of predicted
     trajectories. For nominal and ordinal families this is indexed
     first by boosted subproblem and then by subject.}

  \item{muhat}{Predicted full profiles reconstructed on
    \code{time.grid}, where \code{time.grid} is the fitted training time
    grid used by the model.}
  
  \item{n}{Number of subjects in the prediction data.}

  \item{n.q}{Number of boosted subproblems.}

  \item{ni}{Number of observations or requested time points per subject.}

  \item{nu}{Boosting step size used by the fitted model.}

  \item{nu.vec}{Expanded step-size vector on the time-basis scale.}

  \item{partial}{Logical; whether prediction was requested with a common
  user-supplied time grid.}

  \item{prob.class}{Predicted class probabilities on the original response
  scale for non-continuous families; \code{NULL} for the continuous family.}

  \item{prob.hat.class}{Class probabilities over time for non-continuous
  families; \code{NULL} for the continuous family or when
  \code{use.cv.flag = TRUE}.}

  \item{q.set}{Threshold levels (ordinal) or non-reference levels (binary or
  nominal) defining the boosted subproblems.}

  \item{q.total}{Total number of response levels for non-continuous families.}

  \item{rmse}{Standardized test-set RMSE evaluated at \code{m.opt} when
  prediction responses \code{y} are supplied; otherwise \code{NULL}.}

  \item{time}{A list of observed or requested times for each subject.}

  \item{time.design}{Subject-specific time-design matrices used for
  prediction.}

  \item{time.grid}{The common grid used for prediction; typically the
  training time grid.}

  \item{time.unique}{Sorted unique times appearing in \code{time}.}

  \item{use.cv.flag}{Logical; whether out-of-bag coefficient estimates were
  used.}

  \item{x}{Prediction covariates with one row per subject.}

  \item{x.var.names}{Covariate names expected by the fitted model.}

  \item{y}{Observed prediction-set responses split by subject when
  supplied; otherwise \code{NULL}.}

  \item{y.levels}{Observed response levels from the training fit for
  non-continuous families. \code{NA} for the continuous family.}

  \item{y.mean}{Overall response mean used for standardization.}

  \item{y.org}{Prediction-set responses encoded at the boosted-subproblem
  level when \code{y} is supplied; otherwise \code{NULL}. For continuous and
  binary families this is a subject-level list. For nominal and ordinal
  families it is indexed first by boosted subproblem and then by subject.}

  \item{y.reference}{Reference response level used by the nominal family;
  \code{NULL} otherwise.}

  \item{y.sd}{Overall response standard deviation used for standardization.}
}
}


\author{
Amol Pande, Udaya B. Kogalur and Hemant Ishwaran
}

\references{
Friedman J.H. (2001). Greedy function approximation: a gradient boosting
machine. \emph{The Annals of Statistics}, 29(5):1189--1232.

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
  \code{\link{partial.plot.boostmtree}},
  \code{\link{plot.boostmtree}},
  \code{\link{print.boostmtree}},
  \code{\link{simLong}},
  \code{\link{spirometry}},
  \code{\link{vimp.boostmtree}}
}

\examples{
## -------------------------------------------------------------
## Continuous longitudinal prediction on a held-out test set.
## -------------------------------------------------------------
set.seed(31)
sim.obj <- simLong(n = 20, n.test = 10, n.time = 4, model = 1,
                   family = "continuous")
dta <- sim.obj$data.list
trn <- sim.obj$train.index

fit <- boostmtree(
  x = dta$features[trn, , drop = FALSE],
  tm = dta$time[trn],
  id = dta$id[trn],
  y = dta$y[trn],
  family = "continuous",
  M = 10,
  verbose = FALSE
)

pred.obj <- predict(
  fit,
  x = dta$features[-trn, , drop = FALSE],
  tm = dta$time[-trn],
  id = dta$id[-trn],
  y = dta$y[-trn]
)

print(pred.obj)

## -------------------------------------------------------------
## Predict full profiles for new subjects on the training time grid.
## -------------------------------------------------------------
new.subjects <- dta$features[trn, , drop = FALSE][1:3, ]
pred.obj <- predict(fit, x = new.subjects)
str(pred.obj$mu[[1]], max.level = 1)

## -------------------------------------------------------------
## Predict on a user-supplied common time grid.
## -------------------------------------------------------------
grid.time <- seq(min(dta$time[trn]), max(dta$time[trn]), length.out = 25)
pred.grid <- predict(
  fit,
  x = new.subjects,
  tm = grid.time,
  partial = TRUE
)

str(pred.grid$mu[[1]], max.level = 1)

\donttest{
## -------------------------------------------------------------
## Binary longitudinal prediction.
## -------------------------------------------------------------
set.seed(44)
sim.bin <- simLong(n = 25, n.test = 10, n.time = 4, model = 2,
                   family = "binary")
dta.bin <- sim.bin$data.list
trn.bin <- sim.bin$train.index

fit.bin <- boostmtree(
  x = dta.bin$features[trn.bin, , drop = FALSE],
  tm = dta.bin$time[trn.bin],
  id = dta.bin$id[trn.bin],
  y = dta.bin$y[trn.bin],
  family = "binary",
  M = 10,
  verbose = FALSE
)

pred.bin <- predict(
  fit.bin,
  x = dta.bin$features[-trn.bin, , drop = FALSE],
  tm = dta.bin$time[-trn.bin],
  id = dta.bin$id[-trn.bin],
  y = dta.bin$y[-trn.bin]
)

print(pred.bin)
}
}
\keyword{predict}
\keyword{boosting}
