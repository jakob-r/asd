\name{gtreatsel.sim}
\alias{gtreatsel.sim}
\title{
ASD simulation for treatment selection
}
\description{
Function \code{treatsel.sim} runs simulations for a trial design that tests a number of experimental treatments against a single control treatment group in a seamless adaptive trial. Test treatments are compared to the control treatment using Dunnett's many-to-one testing procedure. An interim analysis is undertaken using an early outcome measure for each treatment (and control). A decision is made on which of the treatments to take forward, using a pre-defined selection rule. Data are simulated for the final outcome measure, and data from the interim and final analyses for the final outcome measure are combined together using either the inverse normal or Fisher combination test, and hypotheses tested at the selected level.
This function should not generally be called by the user. The more user-friendly function \code{\link{treatsel.sim}} covers most common applications.
}
\usage{
gtreatsel.sim(z1=c(0,0,0),z2=c(0,0,0),zearly=c(0,0,0),v1=c(1,1,1,1),
           v2=c(1,1,1,1),vearly=c(1,1,1,1),corr=0,weight=0.5,
           nsim=1000,seed=12345678,select=0,epsilon=1,thresh=1,
           level=0.025,ptest=seq(1:length(z1)),fu=FALSE,
           method="invnorm")
}
\arguments{
  \item{z1}{
Vector of test statistics for the final outcome measure based on stage 1 data
}
  \item{z2}{
Vector of test statistics for the final outcome measure based on stage 2 data
}
  \item{zearly}{
Vector of test statistics for the early outcome measure
}
  \item{v1}{
Vector of variances for the final outcome measure based on stage 1 data; in format control treatment variance followed by the test treatment variances
}
  \item{v2}{
Vector of variances for the final outcome measure based on stage 2 data; format as \code{v1}
}
  \item{vearly}{
Vector of variances for the early outcome measure; format as \code{v1}
}
  \item{corr}{
Vector of correlations between the early and final outcome measures for the control and test treatments; format as \code{v1}
}
  \item{weight}{
Weighting between stages 1 and 2; default is for equal weighting (0.5)
}
  \item{nsim}{
Number of simulations (maximum=10,000,000)
}
  \item{seed}{
Seed number
}
  \item{select}{
Selection rule type; 0 = select all treatments, 1 = select maximum, 2 = select maximum two, 3 = select maximum three, 4 = epsilon rule (select means within epsilon of maximum), 5 = randomly select  a single treatment and 6 = threshold rule (select means greater than or equal to threshold). See \code{\link{select.rule}}
}
  \item{epsilon}{
For \code{select} = 4, set epsilon criterion
}
  \item{thresh}{
For \code{select} = 6, set threshold criterion
}
  \item{level}{
Test level (default=0.025)
}
  \item{ptest}{
Vector of treatment numbers for determining power; for example, c(1,2) will count rejections of one or both hypotheses for testing treatments 1 and 2 against control
}
  \item{fu}{
Logical indicating whether patients from dropped treatments (after interim selection) should be followed-up; default \code{TRUE}
}
  \item{method}{
Select combination method; available options are \dQuote{\code{invnorm}} or \dQuote{\code{fisher}}, with default \dQuote{\code{invnorm}}
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
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{treatsel.sim}}
}
\examples{

gtreatsel.sim(z1=c(1,0,2),z2=c(1,0,2),zearly=c(1,0,1), 
    v1=c(1,1,1,1),v2=c(1,1,1,1),vearly=c(1,1,1,1), 
    corr=0,weight=0.25,nsim=100,seed=12345678,
    select=1,level=0.025,ptest = c(1:3),method="fisher")

}
\keyword{design}
