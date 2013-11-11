\name{gsubpop.sim}
\alias{gsubpop.sim}
\title{
ASD simulation for subpopulation selection
}
\description{
Function \code{subpop.sim}  runs simulations for a trial design that tests an experimental treatment against a single control treatment group in a seamless adaptive trial with co-primary analyses in a pre-defined subgroup and the full population. An interim analysis is undertaken using an early outcome measure and a decision is made on whether to continue with both full and subpopulations, the subpopulation only or the full population, using a pre-defined selection rule. A number of different methods to control the family wise error rate are implemented; (i) the treatment is compared to the control in the subpopulation and full populations using Simes test and the inverse normal combination function used to combine p-values before and after design adaptation, (ii) as (i) but the bivariate normal method of Spiessens and Debois (2010) is used to control the type I error rate, (iii) as (i) but a Bonferroni test is used and (iv) a conditional error function approach using the Spiessens and Debois test. Data are simulated for the early and final outcome measures, subpopulation prevalence and correlation between the final and the early outcomes. 
This function should not generally be called by the user. The more user-friendly function \code{\link{subpop.sim}} covers most common applications.
}
\usage{
gsubpop.sim(z.early=NULL,z1=z1,z2=z2,sprev=sprev,
        corr=NULL,selim=NULL,nsim=nsim,seed=12345678,
        level=level,select="thresh",wt=NULL,method="CT-SD")
}
\arguments{
  \item{z.early}{
Vector of test statistics for early outcome subpopulation and full population i.e. c(sub, full)
}
  \item{z1}{
Vector of test statistics for final outcome subpopulation and full population i.e. c(sub, full)
}
  \item{z2}{
Vector of test statistics for final outcome subpopulation and full population, and subpopulation and full population when both are selected i.e. c(sub only, full only, sub, full)
}
  \item{sprev}{
Subpopulation prevalence
}
  \item{corr}{
Correlation between early and final outcomes
}
  \item{selim}{
Upper and lower limits for the difference between test statistics for the threshold rule
}
  \item{nsim}{
Number of simulations (maximum=10,000,000)
}
  \item{seed}{
Seed number
}
  \item{level}{
Test level (default=0.025)
}
  \item{select}{
Selection rule type; available options are \dQuote{\code{thresh}} and \dQuote{\code{futility}}
}
  \item{wt}{
User set weight for combination test
}
  \item{method}{
Test type; available options are \dQuote{\code{CT-Simes}},  \dQuote{\code{CT-SD}}, \dQuote{\code{CT-Bonferroni}} or \dQuote{\code{CEF}}
}
}
\details{
A structured description of the the methodology and the simulation model is given by Friede \emph{et al.} (2012).
}
\value{
 \item{results}{Table of counts; (i) the number of times the subpopulation, full population or both population are selected (\code{n}), (ii) the number of times the subpopulation is rejected when either it alone or both populations are selected (\code{Hs}), (iii) the number of times the full population is rejected when either it alone or both populations are selected (\code{Hf}), (iv) the number of times both populations are rejected (\code{Hs+Hf}) and (v) the number of times the intersection hypothesis is rejected (\code{Hs+f})}
}
\references{
Spiessens B, Debois M. Adjusted significance levels for subgroup analysis in clinical trials. \emph{Contemporary Clinical Trials} 2010;31:647-656.\cr

Jenkins M, Stone A, Jennison C. An adaptive seamless phase II/III design for oncology trials with subpopulation selection using survival endpoints. \emph{Pharmaceutical Statistics} 2011;10:347-356. \cr

Friede T, Parsons N, Stallard N. A conditional error function approach for subgroup selection in adaptive clinical trials. \emph{Statistics in Medicine} 2012;31:409-4320.\cr

}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{subpop.sim}}
}
\examples{

gsubpop.sim(z.early=c(-1,-1),z1=c(-1,-1),z2=c(-1,0,-1,0),sprev=c(0.5,0.5),
        corr=0.5,selim=c(-0.5,0.5),nsim=1000,seed=12345678,level=0.025,
        select="thresh",wt=0.5,method="CT-SD")

}
\keyword{design}
