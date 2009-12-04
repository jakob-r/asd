\name{asd.sim}
\alias{asd.sim}
\title{
ASD Simulation for Normal Data
}
\description{
Function \code{asd.sim} runs simulations for a trial design that tests a number of experimental treatments against a single control treatment group in a seamless adaptive trial. Test treatments are compared to the control treatment using Dunnett\'s many-to-one testing procedure. An interim analysis is undertaken using an early outcome measure, assumed to be normally distributed, and characterized by standardized treatment effects with variances assumed to be equal to one, for each treatment (and control). A decision is made on which of the treatments to take forward, using a pre-defined selection rule. Data are simulated for the final outcome measure, also characterized by standardized treatment effects, with variance assumed equal to one and a fixed correlation between the final and the early outcomes. Data from the interim and final analyses for the final outcome measure are combined together using either the inverse normal or Fisher combination test and hypotheses are tested at the selected level.
}
\usage{
asd.sim(nsamp = c(32, 32), early = c(0, 0, 0), final = c(0, 0, 0),
            nsim = 1000, corr = 0, seed = 12345678, select = 0, epsilon = 1,
            thresh = 1, level = 0.025, ptest = seq(1:length(early)),
            reall = FALSE, fu = FALSE, estim = FALSE, method = "invnorm")
}
\arguments{
  \item{nsamp}{
Vector of sample sizes for each treatment group, at interim and final analysis or a vector of the total sample sizes at each stage if reallocate is TRUE for the epsilon and futility rule.
}
  \item{early}{
Vector of standardized treatment effects for the early outcome measure (maximum 8 treatments).
}
  \item{final}{
Vector of standardized treatment effects for the final outcome measure.
}
  \item{nsim}{
Number of simulations (maximum=10,000,000).
}
  \item{corr}{
Correlation between early and final outcomes.
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
Vector of treatment numbers for determining power. For example, c(1,2) will count rejections of one or both hypotheses for testing treatments 1 and 2 against the control.
}
  \item{reall}{
Logical indicating if reallocation should be used for the second stage. This is implemented only for the epsilon and threshold rules. Default is FALSE.
}
  \item{fu}{
Logical indicating whether patients from dropped treatments (after interim selection) should be followed-up. Default is FALSE.
}
  \item{estim}{
Logical indicating if mean and variance of final outcome should be estimated. Default is FALSE.
}
  \item{method}{
Select combination method. Available options are \dQuote{\code{invnorm}} or \dQuote{\code{fisher}}. Default is \dQuote{\code{invnorm}}.
}
}
\details{
A structured description of the algorithm used to perform the simulations is given below:\cr
(1) Generate a random sample of bivariate normal effects, using the function \code{\link{simeans.binormal}}, based on a treatment group sample size of \code{nsamp} for the interim and final analyses for the control and \emph{n} test treatments with expected standardized effects given by the \code{early} and \code{final} outcomes.\cr
(2) Calculate \emph{n} test statistics for the \code{early} outcome (based on data from (1)) comparing each test treatment to the control and use the function \code{\link{select.rule}} at the interim analysis (using options \code{type}, \code{epsilon}, \code{thresh}, and \code{reall} as appropriate) to select treatments for the final analysis.\cr
(3) Implement the Dunnett test (\code{\link{dunnett.test}}) for the test statistics from (2) (\emph{dunnett1}), with appropriate modification depending on whether patients are followed-up.\cr
(4) Generate a random sample of normal effects based on a treatment group sample size of \code{nsamp}[2] for the final analysis for the control and \emph{n} test treatments with expected standardized effects given by \code{final}. Calculate \emph{n} test statistics comparing each treatment to the control, for the selected treatments from(2). Implement the Dunnett test (\code{\link{dunnett.test}}) for these test statistics (\emph{dunnett2}), with appropriate modification depending on which treatments are selected.\cr
(5) Outputs \emph{dunnett1} and \emph{dunnett2} are combined using either the \dQuote{\code{fisher}} or the \dQuote{\code{inverse normal}} combination function (\code{\link{combn.test}}).\cr
(6) Hypotheses are tested using the output of the combination function and function \code{\link{hyp.test}} at the selected \code{level}.\cr
(7) For each iteration an indicator records whether each elementary hypothesis (comparing each test treatment to the control) is rejected or accepted. If one or more of the treatments selected using option \code{ptest} is rejected then this is recorded.\cr
(8) Steps (1) to (7) are repeated \code{nsim} times.\cr
(9) The results are summarised by \emph{count.total} which reports the number of times one or more treatments are selected at interim, \emph{select.total} which reports the number of times each test treatment is selected, \emph{reject.total} which gives the number of times each hypothesis is rejected and \emph{sim.reject} which gives the number of times one or more of the treatments selected using \code{ptest} is rejected.\cr
If \code{reall} is TRUE, then \code{nsamp}[1] and \code{nsamp}[2] represent the total number of patients at each stage, and patients are divided equally between treatment groups depending on the number of test treatments and the number of test treatments selected at the interim analysis. If \code{estim} is TRUE then means and variances are estimated and reported for each of the two stages.\cr
}
\value{
 \item{count.total}{Number of times one or more treatments are selected.}
 \item{select.total}{Number of times each test treatment is selected.}
 \item{reject.total}{Number of times each hypothesis is rejected.}
 \item{sim.reject}{Number of times one or more of the treatments selected using \code{ptest} is rejected.}
 \item{est.mean}{Estimates of the mean standardized treatment effects at stages 1 and 2 for each treatment. When \code{estim} is TRUE.}
 \item{est.N}{Estimates of the denominator population for estimating the variance at stages 1 and 2 for each treatment. When \code{estim} is TRUE.}
 \item{est.var}{Estimates of the variances of the standardized treatment effects at stages 1 and 2 for each treatment. When \code{estim} is TRUE.}
}
\references{
Friede T, Parsons N, Stallard S, Todd S, Valdes Marquez E, Chataway J, Nicholas R. Designing a Seamless Phase II/III Clinical Trial using Early Outcomes for Treatment Selection: an Application in Multiple Sclerosis. Submitted to \emph{Statistics in Medicine}.\cr

Bretz F, Schmidli H, Koenig F, Racine A, Maurer W. Confirmatory seamless phase II/III clinical trials with hypotheses selection at interim: General concepts. \emph{Biometrical Journal} 2006;48:623-634.\cr
}
\author{
Nick Parsons (\email{nick.parsons@warwick.ac.uk})
}
\seealso{
\code{\link{gasd.sim}}
}
\examples{

## 3 treatments, one effective with size 0.3, and one effective with size 0.2
## Correlation between early and final outcomes is 0.3
## Select one treatment at the interim analysis
## 32 patients in each group at stage 1 and stage 2
## Count rejections in treatment groups 1 and/or 2
asd.sim(nsamp=c(32,32),early=c(0.3,0.2,0),final=c(0.3,0.2,0),
       nsim=100,corr=0.3,seed=145514,select=1,epsilon=1,reall=FALSE,
       level=0.025,ptest=c(1,2),fu=FALSE,method="invnorm")

## Allow more than one treatment to be selected at the interim analysis
asd.sim(nsamp=c(32,32),early=c(0.3,0.2,0),final=c(0.3,0.2,0),
       nsim=100,corr=0.3,seed=145514,select=4,epsilon=0.5,reall=FALSE,
       level=0.025,ptest=c(1,2),fu=FALSE,method="invnorm")

}
\keyword{design}
