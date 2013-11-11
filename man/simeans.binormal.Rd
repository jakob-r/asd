\name{simeans.binormal}
\alias{simeans.binormal}
\title{
Simulate Bivariate Normal Means
}
\description{
Simulates bivariate normal means; for use with \code{asd.sim} and \code{gasd.sim} in ASD.
}
\usage{
simeans.binormal(n = n, means = means, vars = vars, corr = corr)
}
\arguments{
  \item{n}{
Number of records used to calculate means
}
  \item{means}{
Vector of expected means for two samples
}
  \item{vars}{
Vector of expected variances for two samples
}
  \item{corr}{
Correlation between two samples
}
}
\details{
Uses function \code{rmvnorm} from package \code{mvtnorm} to generate means from correlated normal variates.
}
\value{
 \item{samp1}{Mean of sample 1}
 \item{samp2}{Mean of sample 2}
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{treatsel.sim}}, \code{\link{dunnett.test}}, \code{\link{hyp.test}}, \code{\link{select.rule}}, \code{\link{combn.test}}
}
\examples{

# need to load mvtnorm
library(mvtnorm)

# generate data
set.seed(1234)
simeans.binormal(n=10,means=c(2,3),vars=c(1,5),corr=0.5)

}
\keyword{design}
