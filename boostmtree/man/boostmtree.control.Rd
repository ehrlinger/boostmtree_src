\name{boostmtree.control}
\alias{boostmtree.control}
\title{Create a control object for \code{boostmtree}}
\description{
Construct a control object used by \code{boostmtree()} to manage resampling,
split search, reproducibility, and a few advanced fitting options.
}
\usage{
boostmtree.control(
  ntree = 1,
  bootstrap = "by.root",
  bootstrap.fraction = 0.632,
  sample.matrix = NULL,
  nsplit = NULL,
  samptype = "swor",
  xvar.wt = NULL,
  case.wt = NULL,
  seed = NULL,
  cv.lambda = FALSE,
  cv.rho = TRUE
)
}
\arguments{

  \item{ntree}{Number of trees requested. The current implementation
    supports \code{ntree = 1} only.}
  
  \item{bootstrap}{Resampling rule passed to the underlying tree learner.
    Recognized values are \code{"by.root"}, \code{"by.user"}, and
    \code{"none"}. The default \code{"by.root"} creates out-of-bag subjects and
    is the recommended choice whenever \code{cv.flag = TRUE} or grow-object
    variable importance will be used.}

  \item{bootstrap.fraction}{Sampling fraction used to build the in-bag sample
    when \code{bootstrap = "by.user"} and \code{sample.matrix} is not supplied.}

  \item{sample.matrix}{Optional user-supplied in-bag sampling matrix with
    dimensions \code{n x M}. Entry \code{[i, m]} gives the in-bag count for subject
    \code{i} at boosting iteration \code{m}. When this is used together with
    \code{cv.flag = TRUE}, every column must contain at least one zero so that
    out-of-bag subjects are available.}

  \item{nsplit}{Optional limit on the number of candidate split points examined
    at each node. Smaller values can reduce computation.}

  \item{samptype}{Sampling type passed to the underlying tree learner.}

  \item{xvar.wt}{Optional split-variable weights.}

  \item{case.wt}{Optional case weights.}

  \item{seed}{Optional integer seed used to reproduce the random parts of tree
    fitting. This affects reproducibility only; it does \emph{not} turn
    out-of-bag sampling on or off.}

  \item{cv.lambda}{Logical; should the CV path use the cross-validated mean when
    updating \code{lambda}?}

  \item{cv.rho}{Logical; should the CV path use the cross-validated mean when
    updating \code{rho} and \code{phi}?}

}

\details{

The control object separates routine model arguments such as \code{M},
\code{k}, and \code{nu} from lower-level options that control how the
tree base learner is grown.

The most important option is \code{bootstrap}. In this package, out-of-bag
subjects are used for the cross-validation path stored when
\code{cv.flag = TRUE}, and for grow-object variable importance returned by
\code{vimp.boostmtree()}. The default \code{bootstrap = "by.root"} is therefore
usually the right choice. Setting \code{bootstrap = "none"} disables OOB
sampling. That can be useful for deterministic in-bag fits and for workflows
that rely on a separate held-out test set, but it is not suitable for OOB CV.

When \code{cv.flag = TRUE}, \code{boostmtree()} checks the requested control
object before fitting. If the requested resampling rule would yield no
out-of-bag subjects, the fit is adjusted back to the default OOB rule or the
function stops with a clear message when a user-supplied sampling matrix is not
compatible with OOB CV.

}

\value{
An object of class \code{"boostmtree.control"}.
}

\seealso{
\code{\link{boostmtree}}
}

\examples{
## Default OOB-producing control object.
boostmtree.control(seed = 7)

## Deterministic in-bag fitting with no OOB subjects.
boostmtree.control(bootstrap = "none", seed = 7)
}
