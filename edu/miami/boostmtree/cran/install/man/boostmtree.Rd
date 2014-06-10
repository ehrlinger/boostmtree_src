\name{boostmtree}
\alias{boostmtree}
\alias{_PROJECT_PACKAGE_NAME_}
\title{Boosted multivariate trees for longitudinal data.}
\description{
  A random forest (Breiman, 2001) is grown using user supplied training
  data.  Applies when the response (outcome) is numeric, categorical
  (factor), or right-censored (including competing risk), and yields
  regression, classification, and survival forests, respectively.  The
  resulting forest, informally referred to as a RF-SRC object, contains
  many useful values which can be directly extracted by the user and/or
  parsed using additional functions (see the examples below).  This is
  the main entry point to the \pkg{_PROJECT_PACKAGE_NAME_} package.

}
\usage{
 boostmtree(x,
   tm,
   id,
   y,
   M = 200,
   M.burn = 0,
   nu = 0.01,
   K = 3,
   nknots = 10,
   d = 3,
   lambda,
   lambda.iter = 5,
   pen.ord = 2,
   svd.tol = 1e-6,
   lambda.max = 1e3,
   forest.tol = 1e-3,
   verbose = TRUE,
   ...)
}
\arguments{
    \item{x}{....}
    \item{tm}{....}
    \item{id}{....}
    \item{y}{....}
    \item{M}{...}
    \item{M.burn}{....}
    \item{nu}{....}
    \item{K}{....}
    \item{nknots}{....}
    \item{d}{....}
    \item{lambda}{....}
    \item{lambda.iter}{....}
    \item{pen.ord}{....}
    \item{svd.tol}{....}
    \item{lambda.max}{....}
    \item{forest.tol}{....}
    \item{verbose}{....}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  blah blah blah blah
}
\value{
     An object of class \code{(rfsrc, grow)} with the following
     components:

     \item{call}{The original call to \command{rfsrc}.}
     \item{formula}{The formula used in the call.}
     \item{family}{The family used in the analysis.}
     \item{n}{Sample size of the data (depends upon \code{NA}'s, see \option{na.action}).}
     \item{ntree}{Number of trees grown.}
     \item{mtry}{Number of variables randomly selected for splitting at each node.}
     \item{nodesize}{Minimum size of terminal nodes.}
     \item{nodedepth}{Maximum depth allowed for a tree.}
     \item{splitrule}{Splitting rule used.}
     \item{nsplit}{Number of randomly selected split points.}
     \item{yvar}{y-outcome values.}
     \item{yvar.names}{A character vector of the y-outcome names.}
     \item{xvar}{Data frame of x-variables.}
     \item{xvar.names}{A character vector of the x-variable names.}
     \item{xvar.wt}{Vector of non-negative weights specifying the
        probability used to select a variable for splitting a node.}
     \item{split.wt}{Vector of non-negative weights where entry
        \code{k}, after normalizing, is the multiplier by which the
	  split statistic for a covariate is adjusted.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
         forest. Vector of length \option{ntree}.  A value of zero indicates
         a rejected tree (can occur when imputing missing data).
         Values of one indicate tree stumps.}
     \item{proximity}{Proximity matrix recording the frequency of pairs of data points
         occur within the same terminal node.}
     \item{forest}{If \option{forest=TRUE}, the forest object is returned.
         This object is used for prediction with new test data
	 sets and is required for other R-wrappers.}
     \item{membership}{Matrix recording terminal node membership where
       each column contains the node number that a case falls in for that tree.}
     \item{inbag}{Matrix recording inbag membership where each column
     contains the number of times that a case appears in the bootstrap
     sample for that tree.}
     \item{var.used}{Count of the number of times a variable is used in
       growing the forest.}
     \item{imputed.indv}{Vector of indices for cases with missing
         values.}
     \item{imputed.data}{Data frame of the imputed data. The first
	 column(s) are reserved for the y-responses, after which the
	 x-variables are listed.}
     \item{split.depth}{Matrix [i][j] or array [i][j][k] recording the
           minimal depth for variable [j] for case [i], either averaged over
           the forest, or by tree [k].}
     \item{node.stats}{Split statistics returned when
     \code{statistics=TRUE} which can be parsed using \command{stat.split}.}
     \item{err.rate}{Tree cumulative OOB error rate.}
     \item{importance}{Variable importance (VIMP) for each x-variable.}
     \item{predicted}{In-bag predicted value.}
     \item{predicted.oob}{OOB predicted value.}\cr
}
\note{
  \enumerate{
    \item The returned object depends heavily on the family.  In
    particular, \code{predicted} and \code{predicted.oob} are the
    following values calculated using in-bag and OOB data.
    \enumerate{
      \item For regression, a vector of predicted y-responses.

      \item For classification, a matrix with columns containing the
      estimated class probability for each class.

      \item For survival, a vector of mortality values (Ishwaran et al.,
      2008) representing estimated risk for each individual calibrated
      to the scale of the number of events (as a specific example, if
      \emph{i} has a mortality value of 100, then if all individuals had
      the same x-values as \emph{i}, we would expect an average of 100
      events).  Also included in the grow object are matrices containing
      the CHF and survival function.  Each row corresponds to an
      individual's ensemble CHF or survival function evaluated at each
      time point in \code{time.interest}.

      \item For competing risks, a matrix with one column for each event
      recording the expected number of life years lost due to the event
      specific cause up to the maximum follow up (Ishwaran et al.,
      2013).  The grow object also contains the cause-specific
      cumulative hazard function (CSCHF) and the cumulative incidence
      function (CIF) for each event type.  These are encoded as a
      three-dimensional array, with the third dimension used for the
      event type, each time point in \code{time.interest} making up the
      second dimension (columns), and the case (individual) being the
      first dimension (rows).
    }

    \item Different R-wrappers are provided to aid in parsing the grow
    object.
  }
}
\author{
    Hemant Ishwaran and Udaya B. Kogalur
}
\references{
    Breiman L., Friedman J.H., Olshen R.A. and Stone C.J.
    \emph{Classification and Regression Trees}, Belmont, California, 1984.

    Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    Cutler A. and Zhao G. (2001). Pert-Perfect random tree ensembles.
    \emph{Comp. Sci. Statist.}, 33: 490-497.

    Gray R.J. (1988).  A class of k-sample tests for comparing the
    cumulative incidence of a competing risk, \emph{Ann. Statist.},
    16: 1141-1154.

    Harrell et al. F.E. (1982).  Evaluating the yield of medical tests,
    \emph{J. Amer. Med. Assoc.}, 247:2543-2546.

    Hothorn T. and Lausen B. (2003). On the exact distribution of maximally selected
    rank statistics, \emph{Comp. Statist. Data Anal.}, 43:121-137.

    Ishwaran H. (2007).  Variable importance in binary regression
    trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.

    Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7(2):25-31.

    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.

    Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and
    Lauer M.S. (2010).  High-dimensional variable selection for survival
    data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.

    Ishwaran H. (2013).  The effect of splitting on random forests.

    Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2014). Random survival forests for competing risks.  To appear
    in \emph{Biostatistics}.

    Lin Y. and Jeon Y. (2006).  Random forests and adaptive nearest
    neighbors, \emph{J. Amer. Statist. Assoc.}, 101:578-590.

    LeBlanc M. and Crowley J. (1993).  Survival trees by goodness of split,
    \emph{J. Amer. Statist. Assoc.}, 88:457-467.

    Loh W.-Y and Shih Y.-S (1997).  Split selection methods for
    classification trees, \emph{Statist. Sinica}, 7:815-840.

    Mogensen, U.B, Ishwaran H. and Gerds T.A. (2012). Evaluating random
    forests for survival analysis using prediction error curves,
    \emph{J.  Statist. Software}, 50(11): 1-23.

    Segal M.R. (1988).  Regression trees for censored data,
    \emph{Biometrics},  44:35-47.
}
\seealso{
  \command{\link{find.interaction}},
  \command{\link{impute.rfsrc}},
  \command{\link{max.subtree}},
  \command{\link{plot.competing.risk}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival}},
  \command{\link{plot.variable}},
  \command{\link{predict.rfsrc}},
  \command{\link{print.rfsrc}},
  \command{\link{rf2rfz}},RF_MLTVR_BEG
  \command{\link{rfsrcSyn}},
  \command{\link{stat.split}},RF_MLTVR_END
  \command{\link{var.select}},
  \command{\link{vimp}}
}
\examples{
}
\keyword{forest}

