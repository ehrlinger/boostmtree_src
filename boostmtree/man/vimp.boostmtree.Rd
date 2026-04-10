\name{vimp.boostmtree}
\alias{vimp.boostmtree}
\title{Permutation Variable Importance for Boosted Tree Models}

\description{
Compute permutation-based variable importance for fitted \code{boostmtree}
objects and for prediction objects produced by \code{predict.boostmtree()}.
}

\usage{
vimp.boostmtree(object, x.names = NULL, joint = FALSE)
}

\arguments{

  \item{object}{
    A fitted object of class \code{(boostmtree, grow)} or a prediction object of
    class \code{(boostmtree, predict)}.}

  \item{x.names}{
    Optional character vector naming the covariates to assess. If omitted,
    importance is computed for all available covariates.}

  \item{joint}{
    Logical value indicating whether the variables listed in \code{x.names} should
    be permuted jointly. When \code{FALSE}, each variable is permuted separately.}

}

\details{

Variable importance is computed by permuting one or more covariates and
then measuring how much predictive accuracy deteriorates.

For grow objects, the procedure uses the out-of-bag prediction path
stored in the fitted object. This requires two things: the model must
have been fit with \code{cv.flag = TRUE}, and the resampling rule must
have produced out-of-bag subjects at every boosting iteration used by
the importance calculation. In ordinary use, this happens automatically
because the default control object uses \code{bootstrap = "by.root"}. If
a model was fit with \code{bootstrap = "none"}, then grow-object
variable importance is not available because there are no OOB subjects
to perturb.

For prediction objects, the procedure uses the supplied test responses
and reports the relative increase in test-set RMSE after
permutation. This route does \emph{not} require OOB sampling because the
comparison is made on the held-out prediction data.

In longitudinal settings, the returned object separates three effects:
\describe{
  \item{\code{main}}{importance of the baseline covariate effect}
  \item{\code{interaction}}{importance of the covariate-time interaction}
  \item{\code{time.effect}}{importance of the time basis alone}
}

The returned value has class \code{vimp.boostmtree}. It can be plotted
directly with \code{plot()}.
}
\value{
An object of class \code{vimp.boostmtree}. Its main components are:
\describe{
  \item{\code{main}}{matrix of permutation importance values for the main covariate effects}
  \item{\code{interaction}}{matrix of time-interaction importance values for longitudinal fits, or \code{NULL}}
  \item{\code{time.effect}}{vector containing the importance of the time basis alone, or \code{NULL}}
  \item{\code{x.var.names}}{names of the assessed covariates}
  \item{\code{metric}}{description of the accuracy measure used in the comparison}
}
}

\references{
Friedman, J. H. (2001). Greedy function approximation: A gradient boosting
machine. \emph{Annals of Statistics}, 29, 1189--1232.

Pande, A., Ishwaran, H., Blackstone, E. H., Rajeswaran, J., and Gillinov, M.
(2022). Application of gradient boosting in evaluating surgical ablation for
atrial fibrillation. \emph{SN Computer Science}, 3, 466.

Pande, A., Ishwaran, H., and Blackstone, E. H. (2022). Boosting for
multivariate longitudinal responses. \emph{SN Computer Science}, 3, 186.
}

\examples{
\donttest{
## -------------------------------------------------------------
## Variable importance from a fitted continuous longitudinal model.
## For grow objects this requires cv.flag = TRUE. The default
## control object already provides OOB subjects.
## -------------------------------------------------------------
set.seed(19)
sim.obj <- simLong(n = 100, n.time = 4, model = 2, family = "continuous")
dta <- sim.obj$data.list

fit <- boostmtree(
  x = dta$features,
  tm = dta$time,
  id = dta$id,
  y = dta$y,
  family = "continuous",
  M = 50,
  cv.flag = TRUE,
  verbose = TRUE
)

vimp.obj <- vimp.boostmtree(fit, x.names = c("x1", "x2"))
plot(vimp.obj)


## -------------------------------------------------------------
## Variable importance from a held-out test set.
## This route does not rely on OOB sampling because the comparison
## is made on the prediction object.
## -------------------------------------------------------------
set.seed(23)
sim.obj <- simLong(n = 200, n.test = 100, n.time = 4, model = 2,
                   family = "continuous")
dta <- sim.obj$data.list
trn <- sim.obj$train.index

fit <- boostmtree(
  x = dta$features[trn, , drop = FALSE],
  tm = dta$time[trn],
  id = dta$id[trn],
  y = dta$y[trn],
  family = "continuous",
  M = 50,
  verbose = TRUE,
  control = boostmtree.control(bootstrap = "none", seed = 23)
)

pred.obj <- predict(
  fit,
  x = dta$features[-trn, , drop = FALSE],
  tm = dta$time[-trn],
  id = dta$id[-trn],
  y = dta$y[-trn]
)

vimp.test <- vimp.boostmtree(pred.obj, x.names = c("x1", "x2"))
plot(vimp.test)
}
}
