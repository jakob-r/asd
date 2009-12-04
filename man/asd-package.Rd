\name{asd-package}
\alias{asd-package}
\alias{asd}
\docType{package}
\title{Simulation Tools for Adaptive Seamless Design (ASD)}
\description{
Functions to run simulations for a trial design that tests a number of experimental treatments against a single control treatment group in a seamless adaptive trial. Test treatments are compared to the control treatment using Dunnett\'s many-to-one testing procedure. An interim analysis is undertaken using an early outcome measure. A decision is made on which of the treatments to take forward, using a pre-defined selection rule. Data are simulated for the final outcome measure that is allowed to be correlated with the early outcome measure. Data from the interim and final analyses for the final outcome measure are combined together using either the inverse normal or Fisher combination test and hypotheses are either rejected or accepted after controlling the familywise error rate at the selected level.
}
\details{
\tabular{ll}{
Package: \tab asd\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-12-04\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Simulations are run using the functions \code{\link{asd.sim}} and \code{\link{gasd.sim}}.
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\references{
Adaptive designs are described in more detail:\cr

Thall PF, Simon R, Ans Ellenberg SS. A two-stage design for choosing amongst several experimental treatments and a control in clinical trials. \emph{Biometrics} 1988;45:537-547.\cr

Thall PF, Simon, R, Ans Ellenberg SS. Two-stage selection and testing designs for comparative clinical trials. \emph{Biometrika} 1989;75,303-310.\cr

Bauer P, Kieser M. Combining different phases in the development of medical treatments within a single trial. \emph{Statistics in Medicine} 1999;18:1833-1848.\cr

Stallard N, Todd S. Sequential designs for phase II and phase III clinical trials incorporating treatment selection. \emph{Statistics in Medicine} 2003;22:689-703.\cr

Posch M, Koenig F, Branson M, Brannath W, Dunger-Baldauf C, Bauer P. Testing and estimation in flexible group sequential designs with adaptive treatment selection. \emph{Statistics in Medicine} 2005;24:3697-3714.\cr

Bretz F, Schmidli H, Koenig F, Racine A, Maurer W. Confirmatory seamless phase II/III clinical trials with hypotheses selection at interim: General concepts. \emph{Biometrical Journal} 2006;48:623-634.\cr

Koenig F, Brannath W, Bretz F, Posch M. Adaptive Dunnett tests for treatment selection. \emph{Statistics in Medicine} 2008;27:1612-1625.\cr

Stallard N, Friede T. A group-sequential design for clinical trials with treatment selection. \emph{Statistics in Medicine} 2008;27:6209-6227.\cr

}
\seealso{
\code{\link{asd.sim}}, \code{\link{gasd.sim}}
}
\keyword{design}