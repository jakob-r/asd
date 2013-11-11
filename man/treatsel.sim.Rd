\name{treatsel.sim}
\alias{treatsel.sim}
\title{
ASD simulation for treatment selection
}
\description{
Function \code{treatsel.sim} runs simulations for a trial design that tests a number of experimental treatments against a single control treatment group in a seamless adaptive trial. Test treatments are compared to the control treatment using Dunnett's many-to-one testing procedure. An interim analysis is undertaken using an early outcome measure for each treatment (and control). A decision is made on which of the treatments to take forward, using a pre-defined selection rule. Data are simulated for the final outcome measure, and data from the interim and final analyses for the final outcome measure are combined together using either the inverse normal or Fisher combination test, and hypotheses tested at the selected level.
}
\usage{
treatsel.sim(n=list(stage1=32,stage2=32),
             effect=list(early=c(0,0,0),final=c(0,0,0)),
             outcome=list(early="N",final="N"),nsim=1000,
             corr=0,seed=12345678,select=0,epsilon=1,
             weight=NULL,thresh=1,level=0.025,ptest=c(1),
             method="invnorm",fu=FALSE,file = "")
}
\arguments{

  \item{n}{
List giving sample sizes for each treatment group at stage 1 (interim) and stage 2 (final) analyses
}
  \item{effect}{
List giving effect sizes for early and final outcomes
}
  \item{outcome}{
List giving outcome type for early and final outcomes; available options are \dQuote{\code{N}}, \dQuote{\code{T}} and \dQuote{\code{B}}, for normal, time-to-event and binary data
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
Selection rule type (\code{\link{select.rule}}); 0 = select all treatments, 1 = select maximum, 2 = select maximum two, 3 = select maximum three, 4 = epsilon rule (select means within epsilon of maximum), 5 = randomly select  a single treatment and 6 = threshold rule (select means greater than or equal to threshold)
}
  \item{epsilon}{
For \code{select} = 4, set epsilon criterion
}
  \item{weight}{
Optional user set weight for combination test; default is to use those suggested by Jenkins \emph{et al.} (2011)
}
  \item{thresh}{
For \code{select} = 6, set threshold criterion
}
  \item{level}{
Test level (default=0.025)
}
  \item{ptest}{
Vector of treatment numbers for determining power; for example, c(1,2) will count rejections of one or both hypotheses for testing treatments 1 and 2 against the control
}
  \item{method}{
Select combination method; available options are \dQuote{\code{invnorm}} or \dQuote{\code{fisher}}, with default \dQuote{\code{invnorm}}.
}
  \item{fu}{
Logical indicating whether patients from dropped treatments (after interim selection) should be followed-up; default \code{FALSE}
}
  \item{file}{
File name to dump output; if unset will default to R console
}
}
\details{
A structured description of the the methodology and the simulation model is given by Friede \emph{et al.} (2011) and implementation by Parsons \emph{et al.} (2012).
}
\value{
 \item{count.total}{Number of times one or more treatments are selected}
 \item{select.total}{Number of times each test treatment is selected}
 \item{reject.total}{Number of times each hypothesis is rejected}
 \item{sim.reject}{Number of times one or more of the treatments selected using \code{ptest} is rejected}
}
\references{
Friede T, Parsons N, Stallard N, Todd S, Valdes Marquez E, Chataway J, Nicholas R. Designing a Seamless Phase II/III Clinical Trial using Early Outcomes for Treatment Selection: an Application in Multiple Sclerosis. \emph{Statistics in Medicine} 2011;30:1528-1540.\cr

Parsons N, Friede T, Todd S, Valdes Marquez E, Chataway J, Nicholas R, Stallard N. An R package for implementing simulations for seamless phase II/III clinicals trials using early outcomes for treatment selection. \emph{Computational Statistics and Data Analysis} 2012;56:1150-1160.\cr

Bretz F, Schmidli H, Koenig F, Racine A, Maurer W. Confirmatory seamless phase II/III clinical trials with hypotheses selection at interim: General concepts. \emph{Biometrical Journal} 2006;48:623-634.\cr
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{gtreatsel.sim}}
}
\examples{

# two test treatment groups
# effect size = 0.3 for group 1
# for both early and final normal outcomes
# correlation = 0.3
# select one treatment only at interim
treatsel.sim(n=list(stage1=100,stage2=300),
        effect=list(early=c(0,0.3,0),final=c(0,0.3,0)),
        outcome=list(early="N",final="N"),
        nsim=100,corr=0.3,seed=145514,select=1,
        level=0.025,ptest=c(1,2),fu=FALSE,
        method="invnorm",file="")

# five test treatment groups
# correlation = 0.3
# flexible selection rule, with epsilon = 1
treatsel.sim(n=list(stage1=100,stage2=300),
        effect=list(early=c(0,0.3,0.2,0.1,0.3,0.05),
        final=c(0,0.2,0.3,0.2,0.1,0.5)),
        outcome=list(early="N",final="N"),
        nsim=200,corr=0.3,seed=145514,select=4,epsilon=1,
        level=0.025,ptest=c(1:5),method="invnorm")

}
\keyword{design}
