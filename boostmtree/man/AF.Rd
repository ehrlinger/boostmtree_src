\name{AF}
\docType{data}
\alias{AF}
\title{Atrial Fibrillation Data}
\description{
Atrial Fibrillation (AF) data from a randomized clinical trial examining
the effect of surgical ablation as a treatment for AF in patients with
persistent or long-standing persistent AF who required mitral valve surgery.
Patients were randomized into two groups: mitral valve surgery with ablation
and mitral valve surgery without ablation. Patients assigned to the ablation
group were further randomized to one of two procedures: pulmonary vein
isolation (PVI) or biatrial maze procedure. All patients were followed weekly
over a 12-month period. The primary outcome is the presence or absence of AF,
measured as a binary longitudinal response. The dataset includes 228 patients,
with a total of 7,949 AF measurements and an average of 35 observations per
patient.
}

\format{
  A list containing four elements:
  \enumerate{
    \item  The 84 patient variables (features).
    \item  Time points (time).
    \item  Unique patient identifier (id).
    \item  Presence or absence of AF (y).
  }
}

\references{
Gillinov A. M., Gelijns A.C., Parides M.K., DeRose J.J.Jr.,
Moskowitz~A.J. et al. Surgical ablation of atrial fibrillation during mitral
valve surgery.

\emph{The New England Journal of Medicine} 372(15):1399--1408, 2015.
}

\examples{data(AF, package = "boostmtree")}

\keyword{datasets}
