\name{gasd.sim}
\alias{gasd.sim}
\title{
General ASD Simulation
}
\description{
Function \code{gasd.sim} is a generalization of function \code{\link{asd.sim}} for non-normal data through the use of expected test statistics, that are assumed to be normally distributed, and a more general statement of the variances and correlations. The structure of the outputs of this function are equivalent to \code{asd.sim}, but the inputs differ somewhat.
}
\usage{
gasd.sim(z1 = c(0, 0, 0), z2 = c(0, 0, 0), zearly = c(0, 0, 0), v1 = c(1, 1, 1, 1),
           v2 = c(1, 1, 1, 1), vearly = c(1, 1, 1, 1), corr = c(0, 0, 0, 0), weight = 0.5,
           nsim = 1000, seed = 12345678, select = 0, epsilon = 1, thresh = 1, level = 0.025,
           ptest = seq(1:length(z1)), fu = FALSE, method = "invnorm")
}
\arguments{
  \item{z1}{
Vector of test statistics for the final outcome measure based on stage 1 data (maximum 8 treatments).
}
  \item{z2}{
Vector of test statistics for the final outcome measure based on stage 2 data.
}
  \item{zearly}{
Vector of test statistics for the early outcome measure.
}
  \item{v1}{
Vector of variances for the final outcome measure based on stage 1 data; in format control treatment variance followed by the test treatment variances.
}
  \item{v2}{
Vector of variances for the final outcome measure based on stage 2 data; format as \code{v1}.
}
  \item{vearly}{
Vector of variances for the early outcome measure; format as \code{v1}.
}
  \item{corr}{
Vector of correlations between the early and final outcome measures for the control and test treatments; format as \code{v1}.
}
  \item{weight}{
Weighting between stages 1 and 2; default is for equal weighting (0.5).
}
  \item{nsim}{
Number of simulations (maximum=10,000,000).
}
  \item{seed}{
Seed number.
}
  \item{select}{
Selection rule type; 0 = select all treatments, 1 = select maximum, 2 = select maximum two, 3 = select maximum three, 4 = epsilon rule (select means within epsilon of maximum), 5 = randomly select  a single treatment and 6 = threshold rule (select means greater than or equal to threshold). See \code{\link{select.rule}}.
}
  \item{epsilon}{
For \code{select} = 4, set epsilon criterion.
}
  \item{thresh}{
For \code{select} = 6, set threshold criterion.
}
  \item{level}{
Test level (default=0.025).
}
  \item{ptest}{
Vector of treatment numbers for determining power. For example, c(1,2) will count rejections of one or both hypotheses for testing treatments 1 and 2 against control.
}
  \item{fu}{
Logical indicating whether patients from dropped treatments (after interim selection) should be followed-up. Default is TRUE.
}
  \item{method}{
Select combination method. Available options are \dQuote{\code{invnorm}} or \dQuote{\code{fisher}}. Default is \dQuote{\code{invnorm}}.
}
}
\details{
This function is a generalization of \code{asd.sim}; see \code{\link{asd.sim}} for full details. The algorithm implemented in \code{gasd.sim} is the same as that implemented in \code{asd.sim} except that the full trivariate normal distribution of the test statistics (see Friede \emph{et al.} (2010)) is used to generate samples for treatment selection and hypoptheses testing. For full definitions of the inputs necessary for this function see (see Friede \emph{et al.} (2010)).
}
\value{
 \item{count.total}{Number of times one or more treatments are selected.}
 \item{select.total}{Number of times each test treatment is selected.}
 \item{reject.total}{Number of times each hypothesis is rejected.}
 \item{sim.reject}{Number of times one or more of the treatments selected using \code{ptest} is rejected.}
}
\references{
Friede T, Parsons N, Stallard S, Todd S, Valdes Marquez E, Chataway J, Nicholas R. Designing a Seamless Phase II/III Clinical Trial using Early Outcomes for Treatment Selection: an Application in Multiple Sclerosis. Submitted to \emph{Statistics in Medicine}.\cr

Bretz F, Schmidli H, Koenig F, Racine A, Maurer W. Confirmatory seamless phase II/III clinical trials with hypotheses selection at interim: General concepts. \emph{Biometrical Journal} 2006;48:623-634.\cr
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{asd.sim}}
}
\examples{

## General use of function
gasd.sim(z1=c(1.6,0,0),z2=c(1.6,0,0),zearly=c(0.8,0.8,0),
            v1=c(2,5,1,1),v2=c(2,5,1,1),vearly=c(2,10,1,1),
            corr=c(0.5,0.5,0,0),weight=0.5,nsim=100,
            seed=145514,select=1,epsilon=1,thresh=1,
            level=0.025,ptest=c(1),fu=FALSE,method="invnorm")

## asd.sim
## increase nsim to 10000 for more precise estimate
asd.sim(nsamp=c(32,128),early=c(0.3,0,0),final=c(0.3,0,0),
     nsim=100,corr=0,seed=145514,select=1,epsilon=1,reall=FALSE,
     level=0.025,ptest=c(1,2,3),fu=FALSE,method="invnorm")

## gasd.sim should approximtaely replicate asd.sim
## increase nsim to 10000 for more precise estimate
gasd.sim(z1=c(1.2,0,0),z2=c(2.4,0,0),zearly=c(1.2,0,0),
        v1=c(1,1,1,1),v2=c(1,1,1,1),vearly=c(1,1,1,1),
        corr=c(0,0,0,0),weight=0.2,nsim=100,
        seed=145514,select=1,epsilon=1,thresh=1,
        level=0.025,ptest=c(1,2,3),fu=FALSE,method="invnorm")

}
\keyword{design}
