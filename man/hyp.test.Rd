\name{hyp.test}
\alias{hyp.test}
\title{
Closed Testing for ASD
}
\description{
Implements the closure principle (Marcus \emph{et al.}, 1976) for controlling the familywise type I error rate in ASD.
}
\usage{
hyp.test(comb.test, level = level, full.hyp = FALSE)
}
\arguments{
  \item{comb.test}{
Output from function \code{combn.test}
}
  \item{level}{
Test level (default=0.025)
}
  \item{full.hyp}{
Logical indicating whether the full set of intersection hypotheses should be reported; default \code{FALSE}
}
}
\details{
In order to control the familywise type I error rate in the strong sense at the pre-specified level \eqn{\alpha} the closure principle (Marcus \emph{et al.}, 1976) is applied. This means that an individual null hypothesis is rejected if and only if all intersection hypotheses are also rejected at level \eqn{\alpha}.
}
\value{
 \item{reject}{Matrix indicating whether elementary hypotheses have been rejected}
 \item{all.rejects}{Matrix indicating rejections for each intersection hypothesis, if \code{full.hyp}=\code{TRUE}}
 \item{all.hyp}{Matrix labelling each intersection hypothesis, if \code{full.hyp}=\code{TRUE}}
}
\references{
Marcus R, Peritz E, Gabriel KR. On closed testing procedures with special reference to ordered analysis of variance. \emph{Biometrika} 1976;63:655-660.\cr

}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{treatsel.sim}},\code{\link{dunnett.test}}, \code{\link{combn.test}}, \code{\link{select.rule}}, \code{\link{simeans.binormal}}
}
\examples{

stage1 <- dunnett.test(c(0.75,1.5,2.25))
stage2 <- dunnett.test(c(0.15,1.75,2.15))
comb.test <- combn.test(stage1,stage2,weight=0.5)
hyp.test(comb.test,level=0.025,full.hyp=FALSE)

# more output
hyp.test(comb.test,level=0.025,full.hyp=TRUE)

}
\keyword{design}

