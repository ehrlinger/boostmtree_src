\name{predict.boostmtree}
\alias{predict.boostmtree}
\title{Prediction for Boosted multivariate trees for longitudinal data.}
\description{
   Obtain predicted values using a forest.  Also returns performance
   values if the test data contains y-outcomes.
}
\usage{\method{predict}{boostmtree}(object,
 x,
 tm,
 id,
 y,
 M,
 importance = TRUE,
 proximity = FALSE,
 verbose = TRUE,
 eps = 1e-4,
 ...)
}
\arguments{
  \item{object}{....}
  \item{x}{....}
  \item{tm}{....}
  \item{id}{....}
  \item{y}{....}
  \item{M}{....}
  \item{importance}{....}
  \item{proximity}{....}
  \item{verbose}{....}
  \item{eps}{....}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Predicted values are obtained by dropping test data down the grow
  forest (the forest grown using the training data).  The overall error
  rate and VIMP are also returned if the test data contains y-outcome
  values.  Single as well as joint VIMP measures can be requested.  Note
  that calculating VIMP can be computationally expensive (especially
  when the dimension is high), thus if VIMP is not needed, computational
  times can be significantly improved by setting
  \option{importance="none"} which turns VIMP off.

  Setting \option{na.action="na.impute"} imputes missing
  test data (x-variables and/or y-outcomes).  Imputation uses the
  grow-forest such that only training data is used when imputing test
  data to avoid biasing error rates and VIMP (Ishwaran et al. 2008).
  The \option{na.action="na.random"} implements a cruder, but faster
  missing data algorithm.  See the \command{rfsrc} help file for
  details.

  If no test data is provided, then the original training data is used
  and the code reverts to restore mode allowing the user to restore the
  original grow forest.  This is useful because it gives the user the
  ability to extract outputs from the forest that were not asked for in
  the original grow call.  See the examples below for illustration.

  If \option{outcome="test"}, the predictor is calculated by using
  y-outcomes from the test data (outcome information must be present).
  In this case, the terminal nodes from the grow-forest are recalculated
  using the y-outcomes from the test set.  This yields a modified
  predictor in which the topology of the forest is based solely on the
  training data, but where the predicted value is based on the test
  data.  Error rates and VIMP are calculated by bootstrapping the test
  data and using out-of-bagging to ensure unbiased estimates.  See the
  examples for illustration.
}
\value{
     An object of class \code{(rfsrc, predict)}, which is a list with the
     following components:
     \item{call}{The original grow call to \command{rfsrc}.}
     \item{family}{The family used in the analysis.}
     \item{n}{Sample size of test data (depends upon \code{NA} values).}
     \item{ntree}{Number of trees in the grow forest.}
     \item{yvar}{Test set y-outcomes or original grow y-outcomes if none.}
     \item{yvar.names}{A character vector of the y-outcome names.}
     \item{xvar}{Data frame of test set x-variables.}
     \item{xvar.names}{A character vector of the x-variable names.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
       grow forest.  Vector of length \code{ntree}.}
     \item{proximity}{Symmetric proximity matrix of the test data.}
     \item{forest}{The grow forest.}
     \item{membership}{Matrix recording terminal node membership for the
       test data where each column contains the node number that a
       case falls in for that tree.}
     \item{inbag}{Matrix recording inbag membership for the test data
       where each column contains the number of times that a case
       appears in the bootstrap sample for that tree.}
     \item{var.used}{Count of the number of times a variable was used in
       growing the forest.}
     \item{imputed.indv}{Vector of indices of records in test data with
         missing values.}
     \item{imputed.data}{Data frame comprising imputed test data.  First
       columns are the y-outcomes.}
     \item{split.depth}{Matrix [i][j] or array [i][j][k] recording the
     minimal depth for variable [j] for case [i], either averaged over
     the forest, or by tree [k].}
     \item{err.rate}{Cumulative OOB error rate for the test data if
       y-outcomes are present.}
     \item{importance}{Test set variable importance (VIMP).  Can be
       \code{NULL}.}
     \item{predicted}{Test set predicted value.}
     \item{predicted.oob}{OOB predicted value (\code{NULL} unless
       \option{outcome="test"}).}\cr

\item{...... class}{for classification settings, additionally the following ......} \cr
      \item{class}{In-bag predicted class labels.}
      \item{class.oob}{OOB predicted class labels (\code{NULL} unless \option{outcome="test"}).}\cr

\item{...... surv}{for survival settings, additionally the following ......} \cr

     \item{chf}{Cumulative hazard function (CHF).}
     \item{chf.oob}{OOB CHF (\code{NULL} unless \option{outcome="test"}).}
     \item{survival}{Survival function.}
     \item{survival.oob}{OOB survival function (\code{NULL} unless \option{outcome="test"}).}
     \item{time.interest}{Ordered unique death times.}
     \item{ndead}{Number of deaths.}\cr

\item{...... surv-CR}{for competing risks, additionally the following ......} \cr

     \item{chf}{Cause-specific cumulative hazard function (CSCHF)
     for each event.}
     \item{chf.oob}{OOB CSCHF for each event (\code{NULL} unless \option{outcome="test"}).}
     \item{cif}{Cumulative incidence function (CIF) for each event.}
     \item{cif.oob}{OOB CIF for each event (\code{NULL} unless \option{outcome="test"}).}
     \item{time.interest}{Ordered unique event times.}
     \item{ndead}{Number of events.}
}
\note{
  The dimensions and values of returned objects depend heavily on the
  underlying family and whether y-outcomes are present in the test data.
  In particular, items related to performance will be \code{NULL} when
  y-outcomes are not present.

  For detailed definitions of returned values (such as \code{predicted})
  see the \command{rfsrc} help file.
}
\author{
    Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.

  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
}
\seealso{
  \command{\link{plot.competing.risk}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival}},
  \command{\link{plot.variable}},
  \command{\link{vimp}}
}
\examples{

}
\keyword{predict}
\keyword{forest}

