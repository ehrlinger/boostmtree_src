\name{boostmtree-package}
\alias{boostmtree-package}
\docType{package}
\title{
  Boosted multivariate trees for longitudinal data.
}
\description{

  Multivariate boosted tree models for repeated-measures outcomes. The
  method combines gradient boosting, subject-level tree partitioning, and
  smooth time effects to estimate subject-specific mean trajectories or
  class probabilities over time.  Supports continuous, binary, nominal,
  and ordinal longitudinal data, as well as univariate responses.

}

\section{Package Overview}{

  This package contains many useful functions and users should read the
  help file in its entirety for details.  However, we briefly mention
  several key functions that may make it easier to navigate and
  understand the layout of the package.

  \enumerate{
    \item \command{\link{boostmtree}}

    This is the main entry point to the package.  It grows a
    multivariate tree using user supplied training data.  Trees are
    grown using the \pkg{randomForestSRC} R-package.

    \item \command{\link{predict.boostmtree}} (\command{predict})

    Used for prediction.  Predicted values are obtained by dropping the
    user supplied test data down the grow forest.  The resulting object
    has class \code{(rfsrc, predict)}.

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
  \code{\link{boostmtree}},
  \code{\link{boostmtree.control}},
  \code{\link{plot.boostmtree}},
  \code{\link{predict.boostmtree}},
  \code{\link{print.boostmtree}},
  \code{\link{simLong}},
  \code{\link{vimp.boostmtree}},
}

\keyword{package}


