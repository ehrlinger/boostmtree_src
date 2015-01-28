\name{_PROJECT_PACKAGE_NAME_-package}
\alias{_PROJECT_PACKAGE_NAME_-package}
\docType{package}
\title{
  Boosted multivariate trees for longitudinal data.
}
\description{
  This package blah blah blah
  blah blah blah
  blah blah blah
  blah blah blah
}
\section{blah blah blah }{
 This package blah blah blah
}
\section{blah blah blah }{
 This package blah blah blah
}
\section{Package Overview}{

  This package contains many useful functions and users should read the
  help file in its entirety for details.  However, we briefly mention
  several key functions that may make it easier to navigate and
  understand the layout of the package.

  \enumerate{
    \item \command{\link{boostmtree}}

    This is the main entry point to the package.  It grows a
    multivariate tree
    using user supplied training data.  We refer to the resulting object
    as a RF-SRC grow object.  Formally, the resulting object has class
    \code{(rfsrc, grow)}.

    \item \command{\link{predict.boostmtree}} (\command{predict})

    Used for prediction.  Predicted values are obtained by dropping the
    user supplied test data down the grow forest.  The resulting object
    has class \code{(rfsrc, predict)}.

  }
}
\author{
    Hemant Ishwaran and Udaya B. Kogalur
}
\references{
    Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7(2):25-31.

    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.

    Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and
    Lauer M.S. (2010).  High-dimensional variable selection for survival
    data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.

    Ishwaran H., Gerds T.A., Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2014). Random survival forests for competing risks.  To appear
    in \emph{Biostatistics}.
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
  \command{\link{var.select}},
  \command{\link{vimp}}
}
\keyword{package}


