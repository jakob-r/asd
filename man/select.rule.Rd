\name{select.rule}
\alias{select.rule}
\title{
Selection Rules for Interim Analysis in ASD
}
\description{
Function \code{select.rule} provides a number of options for selecting treatments at an interim analysis in ASD.
}
\usage{
select.rule(x, type = 0, epsilon = 1, thresh = 1)
}
\arguments{
  \item{x}{
Vector of test statistics.
}
  \item{type}{
Decision rule type; 0, 1, 2, 3, 4, 5 or 6. See below for details. Default is 0.
}
  \item{epsilon}{
For \code{type} = 4, set epsilon criterion.
}
  \item{thresh}{
For \code{type} = 6, set threshold criterion.
}
}
\details{
There are seven types of selction rule available:\cr
(0) Select all treatments.\cr
(1) Select one treatment; largest value of \code{x}.\cr
(2) Select two treatments; two largest values of \code{x}.\cr
(3) Select three treatments; three largest values of \code{x}.\cr
(4) Epsilon rule; select all \code{x} within \code{epsilon} of maximum.\cr
(5) Randomly select one treatment.\cr
(6) Threshold rule; select all \code{x} larger than \code{thresh}.\cr
}
\value{
 \item{select}{Indicator vector that shows treatments selected (1) or not selected (0).}
 \item{z}{Vector of same length as \code{select} set to -Inf if not selected and 0 otherwise. For use with function \code{dunnett.test}.}
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{asd.sim}}, \code{\link{gasd.sim}}, \code{\link{dunnett.test}}, \code{\link{hyp.test}}, \code{\link{combn.test}}, \code{\link{simeans.binormal}}
}
\examples{

# select maximum treatment
select.rule(x=c(5.3,5.2,1.3,4.5,-1.3),type=4,epsilon=1)

}
\keyword{design}

