\name{subpop.sim}
\alias{subpop.sim}
\title{
ASD simulation for subpopulation selection
}
\description{
Function \code{subpop.sim}  runs simulations for a trial design that tests an experimental treatment against a single control treatment group in a seamless adaptive trial with co-primary analyses in a pre-defined subgroup and the full population. An interim analysis is undertaken using an early outcome measure and a decision is made on whether to continue with both full and subpopulations, the subpopulation only or the full population, using a pre-defined selection rule. A number of different methods to control the family wise error rate are implemented; (i) the treatment is compared to the control in the subpopulation and full populations using Simes test and the inverse normal combination function used to combine p-values before and after design adaptation, (ii) as (i) but the bivariate normal method of Spiessens and Debois (2010) is used to control the type I error rate, (iii) as (i) but a Bonferroni test is used and (iv) a conditional error function approach using the Spiessens and Debois test. Data are simulated for the early and final outcome measures, subpopulation prevalence and correlation between the final and the early outcomes. 
}
\usage{
subpop.sim(n=list(stage1=32,enrich=NULL,stage2=32),
          effect=list(early=c(0,0),final=c(0, 0)),
          outcome=list(early="N",final="N"),
          control=list(early=NULL,final=NULL),sprev=0.5,
          nsim=1000,corr=0,seed=12345678,select="thresh",
          weight=NULL,selim=NULL,level=0.025,method="CT-SD",
          sprev.fixed=TRUE,file="")
}
\arguments{
  \item{n}{
List giving sample sizes for each treatment group at stage 1 (interim) and stage 2 (final) analyses; \code{enrich} allows for sample size modifications if the subgroup only is selected at stage 1
}
  \item{effect}{
List giving effect sizes for early and final outcomes
}
  \item{outcome}{
List giving outcome type for early and final outcomes; available options are \dQuote{\code{N}}, \dQuote{\code{T}} and \dQuote{\code{B}}, for normal, time-to-event and binary data
}
  \item{control}{
Optional list giving effect sizes for early and final outcomes
}
  \item{sprev}{
Subpopulation prevalence
}
  \item{nsim}{
Number of simulations (maximum=10,000,000)
}
  \item{corr}{
Correlation between early and final outcomes
}
  \item{seed}{
Seed number
}
  \item{select}{
Selection rule type;  available options are \dQuote{\code{thresh}} and \dQuote{\code{futility}}
}
  \item{weight}{
Optional user set weight for combination test; default is to use those suggested by Jenkins \emph{et al.} (2011)
}
  \item{selim}{
Upper and lower limits for the difference between test statistics for the threshold rule
}
  \item{level}{
Test level (default=0.025)
}
  \item{method}{
Test type; available options are \dQuote{\code{CT-Simes}},  \dQuote{\code{CT-SD}}, \dQuote{\code{CT-Bonferroni}} or \dQuote{\code{CEF}}
}
  \item{sprev.fixed}{
Logical indicating whether subpopulation prevalence is fixed at each simulation; default \code{TRUE}
}
  \item{file}{
File name to dump output; if unset will default to R console
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
\code{\link{gsubpop.sim}}
}
\examples{


# hazard ratio in subgroup = 0.6 and full population = 0.9
# for both early and final time-to-event outcomes
# subgroup prevalence = 0.3 and correlation = 0.5
# futility stopping rule, with limits 0 and 0
subpop.sim(n=list(stage1=100,enrich=200,stage2=300),
           effect=list(early=c(0.6,0.9),final=c(0.6,0.9)),
           sprev=0.3,outcome=list(early="T",final="T"),nsim=100,
           corr=0.5,seed=1234,select="futility",weight=NULL,
           selim=c(0,0),level=0.025,method="CT-SD",file="")

}
\keyword{design}
