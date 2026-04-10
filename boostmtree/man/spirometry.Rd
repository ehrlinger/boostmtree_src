\name{spirometry}
\docType{data}
\alias{spirometry}
\title{Spirometry Data}
\description{
Lung transplant data from 509 patients who underwent lung transplantation
(LTx) at the Cleveland Clinic, comprising 9,471 longitudinal measurements
of forced expiratory volume in one second (FEV1, expressed as a percentage
of predicted). Twenty-three patient and procedure variables were recorded at
the time of transplant. The primary objectives are to evaluate the temporal
trend of FEV1 following LTx, to identify factors associated with post-LTx
FEV1, and to assess differences in trajectories between single and double
LTx procedures.
}

\format{
  A list containing four elements:
  \enumerate{
    \item  The 23 patient variables (features).
    \item  Time points (time).
    \item  Unique patient identifier (id).
    \item  FEV1-outcomes (y).
  }
}
\references{
  Mason D.P., Rajeswaran J., Li L., Murthy S.C., Su J.W., Pettersson G.B.,
  Blackstone E.H. Effect of changes in postoperative spirometry on
  survival after lung transplantation. \emph{J. Thorac. Cardiovasc. Surg.},
  144:197-203, 2012.
}

\examples{data(spirometry, package = "boostmtree")}

\keyword{datasets}
