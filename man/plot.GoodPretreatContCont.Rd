\name{plot GoodPretreatContCont}
\alias{plot GoodPretreatContCont}
\alias{plot.GoodPretreatContCont}
\title{
Graphically illustrates the theoretical plausibility of finding a good pretreatment predictor in the continuous-continuous case}
\description{
This function provides a plot that displays the frequencies, percentages, or cumulative percentages of \eqn{\rho_{min}^{2}} for a fixed value of \eqn{\delta} (given the observed variances of the true endpoint in the control and experimental treatment conditions and a specified grid of values for the unidentified parameter \eqn{\rho(_{T_{0},T_{1}})}; see \code{\link{GoodPretreatContCont}}). For details, see the online appendix of Alonso et al., submitted.}
\usage{
\method{plot}{GoodPretreatContCont}(x, main, col, Type="Percent", Labels=FALSE, 
Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), \dots)
}
\arguments{
 \item{x}{An object of class \code{GoodPretreatContCont}. See \code{\link{GoodPretreatContCont}}.}
\item{main}{The title of the plot.}
\item{col}{The color of the bins.}
 \item{Type}{The type of plot that is produced. When \code{Type=Freq} or \code{Type=Percent}, the Y-axis shows frequencies or percentages of \eqn{\rho_{min}^{2}}. When \code{Type=CumPerc}, the Y-axis shows cumulative percentages of \eqn{\rho_{min}^{2}}.  Default "Percent".}  
 \item{Labels}{Logical. When \code{Labels=TRUE}, the percentage of \eqn{\rho_{min}^{2}} values that are equal to or larger than the midpoint value of each of the bins are displayed (on top of each bin). Only applies when \code{Type=Freq} or \code{Type=Percent}. Default \code{FALSE}.}
   \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{hist()}.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (submitted). Validating predictors of therapeutic success: a causal inference approach. 
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{GoodPretreatContCont}}
\examples{

# compute rho^2_min in the setting where the variances of T in the control
# and experimental treatments equal 100 and 120, delta is fixed at 50,
# and the grid G={0, .01, ..., 1} is considered for the counterfactual 
# correlation rho_T0T1:

MinPred <- GoodPretreatContCont(T0T0 = 100, T1T1 = 120, Delta = 50,
T0T1 = seq(0, 1, by = 0.01))

# Plot the results (use percentages on Y-axis)
plot(MinPred, Type="Percent")

# Same plot, but add the percentages of ICA values that are equal to or 
# larger than the midpoint values of the bins
plot(MinPred, Labels=TRUE)
}

\keyword{Plausibility of a good pretreatment predictor}
\keyword{Continuous-continuous setting}