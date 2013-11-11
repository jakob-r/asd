\name{asd-package}
\alias{asd-package}
\alias{asd}
\docType{package}
\title{Simulation Tools for Adaptive Seamless Design (ASD)}
\description{
Functions to run simulations for trial designs that either (i) test a number of experimental treatments against a single control treatment group in a seamless adaptive trial or (ii) test an experimental treatment against a single control treatment group in a seamless adaptive trial with co-primary analyses in a pre-defined subgroup and the full population.\cr
In setting (i) test treatments are compared to the control treatment using Dunnett's many-to-one testing procedure, with an interim analysis undertaken using an early outcome measure. A decision is made on which of the treatments to take forward using a pre-defined selection rule. Data are simulated for the final outcome measure that is correlated with the early outcome measure. Data from the interim and final analyses for the final outcome measure are combined together using either the inverse normal or Fisher combination test and hypotheses are either rejected or accepted after controlling the familywise error rate at the selected level.\cr
In setting (ii) an interim analysis is undertaken using an early outcome measure and a decision is made on whether to continue with both full and subpopulations, the subpopulation only or the full population, using a pre-defined selection rule. A number of different methods to control the family wise error rate are implemented. Data are simulated for the early and final outcome measures, subpopulation prevalence and correlation between the final and the early outcomes.
}
\details{
\tabular{ll}{
Package: \tab asd\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2013-11-11\cr
License: \tab GPL-3\cr
}
Simulations are run using the functions (i) \code{\link{treatsel.sim}} and (ii) \code{\link{subpop.sim}}. The other functions are not generally to be called by the user.
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\references{
Some useful references to adaptive designs and more specifically to the methodology described here:\cr

Thall PF, Simon R, Ans Ellenberg SS. A two-stage design for choosing amongst several experimental treatments and a control in clinical trials. \emph{Biometrics} 1988;45:537-547.\cr

Thall PF, Simon, R, Ans Ellenberg SS. Two-stage selection and testing designs for comparative clinical trials. \emph{Biometrika} 1989;75,303-310.\cr

Bauer P, Kieser M. Combining different phases in the development of medical treatments within a single trial. \emph{Statistics in Medicine} 1999;18:1833-1848.\cr

Stallard N, Todd S. Sequential designs for phase II and phase III clinical trials incorporating treatment selection. \emph{Statistics in Medicine} 2003;22:689-703.\cr

Posch M, Koenig F, Branson M, Brannath W, Dunger-Baldauf C, Bauer P. Testing and estimation in flexible group sequential designs with adaptive treatment selection. \emph{Statistics in Medicine} 2005;24:3697-3714.\cr

Bretz F, Schmidli H, Koenig F, Racine A, Maurer W. Confirmatory seamless phase II/III clinical trials with hypotheses selection at interim: General concepts. \emph{Biometrical Journal} 2006;48:623-634.\cr

Koenig F, Brannath W, Bretz F, Posch M. Adaptive Dunnett tests for treatment selection. \emph{Statistics in Medicine} 2008;27:1612-1625.\cr

Stallard N, Friede T. A group-sequential design for clinical trials with treatment selection. \emph{Statistics in Medicine} 2008;27:6209-6227.\cr

Friede T, Parsons N, Stallard N, Todd S, Valdes Marquez E, Chataway J, Nicholas R. Designing a Seamless Phase II/III Clinical Trial using Early Outcomes for Treatment Selection: an Application in Multiple Sclerosis. \emph{Statistics in Medicine} 2011;30:1528-1540.\cr

Parsons N, Friede T, Todd S, Valdes Marquez E, Chataway J, Nicholas R, Stallard N. An R package for implementing simulations for seamless phase II/III clinicals trials using early outcomes for treatment selection. \emph{Computational Statistics and Data Analysis} 2012;56:1150-1160.\cr

Friede T, Parsons N, Stallard N. A conditional error function approach for subgroup selection in adaptive clinical trials. \emph{Statistics in Medicine} 2012;31:4309-4320.\cr

}
\seealso{
\code{\link{treatsel.sim}}, \code{\link{subpop.sim}}
}
\keyword{design}