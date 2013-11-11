\name{dunnett.test}
\alias{dunnett.test}
\title{
Dunnett Test
}
\description{
Implements Dunnett's test (Dunnett, 1955) for many-to-one comparisons.
}
\usage{
dunnett.test(Z = Z, select = rep(1, length(Z)))
}
\arguments{
  \item{Z}{
A vector of test statistics
}
  \item{select}{
A vector of length \code{Z}; to include treatments set values to one and to exclude treatments set values to zero
}
}
\details{
A many-to-one comparison test for the the null hypothesis that all the treatment effects are equal to zero against the alternative that at least one is larger than zero.
}
\value{
 \item{pvalues}{A list of matrices of \emph{p}-values for all intersection hypotheses}
 \item{zscores}{A list of matrices of z-scores for all intersection hypotheses}
 \item{hyp.comb}{A list of matrices indicating the structure of the intersection hypotheses}
}
\references{
Dunnett CW. A multiple comparison procedure for comparing several treatments with a control. \emph{Journal of the American Statistical Association} 1955;50:1096-1121.\cr

}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{treatsel.sim}}, \code{\link{combn.test}}, \code{\link{hyp.test}}, \code{\link{select.rule}}, \code{\link{simeans.binormal}}
}
\examples{

dunnett.test(c(0.75,1.5,2.25))

# select two treatments only
dunnett.test(c(0.75,1.5,2.25),select=c(1,1,0))

# set test statistic to -Inf
dunnett.test(c(0.75,1.5,-Inf))

}
\keyword{design}
