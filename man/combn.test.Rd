\name{combn.test}
\alias{combn.test}
\title{
Combination Tests for ASD
}
\description{
Implements weighted inverse normal and Fisher combination tests for combining \emph{p}-values for adaptive seamless designs.
}
\usage{
combn.test(stage1, stage2, weight = 0.5, method = "invnorm")
}
\arguments{
  \item{stage1}{
Output from function \code{dunnett.test} from stage 1 of an ASD.
}
  \item{stage2}{
Output from function \code{dunnett.test} from stage 2 of an ASD.
}
  \item{weight}{
Weight indicating how \emph{p}-values from stages 1 and 2 are combined. Default weight is 0.5 indicating equal weighting between stages (0<\code{weight}<1).
}
  \item{method}{
Select combination test method. Available options are \dQuote{\code{invnorm}} or \dQuote{\code{fisher}}. Default is \dQuote{\code{invnorm}}.
}
}
\details{
The basic ideas of the combination test approach were proposed by Bauer and Kieser (1999) and make use of a combination function (Bauer and Kohne, 1994) to combine stagewise \emph{p}-values to allow for interim adaptations and the application of the closed test principle (Marcus \emph{et al.}, 1976) to control the overall test size across multiple hypotheses.
}
\value{
 \item{method}{Selected method of combining \emph{p}-values.}
 \item{zscores}{Z-scores for each hypothesis.}
 \item{hyp.comb}{A list of matrices indicating the structure of the intersection hypotheses.}
 \item{weights}{Weights used for each stage.}
}
\references{
Bauer P, Kieser M. Combining different phases in the development of medical treatments within a single trial. \emph{Statistics in Medicine} 1999;18:1833-1848.\cr

Bauer P, Kohne K. Evaluation of experiments with adaptive interim analyses. \emph{Biometrics} 1994;50:1029-1041.\cr

Marcus R, Peritz E, Gabriel KR. On closed testing procedures with special reference to ordered analysis of variance. \emph{Biometrika} 1976;63:655-660.\cr

Lehmacher W, Wassmer G. Adaptive sample size calculations in group sequential trials. \emph{Biometrics} 1999;55:1286-1290.
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{asd.sim}}, \code{\link{gasd.sim}}, \code{\link{dunnett.test}}, \code{\link{hyp.test}}, \code{\link{select.rule}}, \code{\link{simeans.binormal}}
}
\examples{

stage1 <- dunnett.test(c(0.75,1.5,2.25))
stage2 <- dunnett.test(c(0.15,1.75,2.15))
combn.test(stage1,stage2,weight=0.5,method="invnorm")

}
\keyword{design}

