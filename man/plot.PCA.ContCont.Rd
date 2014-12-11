\name{plot PCA.ContCont}
\alias{plot PCA.ContCont}
\alias{plot.PCA.ContCont}
\title{
Plots the Predictive Causal Association in the continuous-continuous case
}
\description{
This function provides a plot that displays the frequencies, percentages, or cumulative percentages of the Predictive Causal Association (PCA; \eqn{\rho_{\psi}}). These figures are useful to examine the sensitivity of the obtained results with respect to the assumptions regarding the correlations between the counterfactuals (for details, see Alonso et al., submitted). Optionally, it is also possible to obtain plots that are useful in the examination of the plausibility of finding a good pretreatment predictor.}
\usage{
\method{plot}{PCA.ContCont}(x, Xlab.PCA, Main.PCA, Type="Percent", 
Labels=FALSE, PCA=TRUE, Good.Pretreat=FALSE, EffectT0T1=FALSE, 
Main.Good.Pretreat, Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), 
col, \dots)
}
\arguments{
  \item{x}{An object of class \code{PCA.ContCont}. See \code{\link{PCA.ContCont}}.}
\item{Xlab.PCA}{The legend of the X-axis of the PCA plot. Default "\eqn{\rho_{\psi}}".}
 \item{Main.PCA}{The title of the PCA plot. Default "PCA".}  
 \item{Type}{The type of plot that is produced. When \code{Type=Freq} or \code{Type=Percent}, the Y-axis shows frequencies or percentages of \eqn{\rho_{\psi}} and/or \eqn{\delta}. When \code{Type=CumPerc}, the Y-axis shows cumulative percentages of \eqn{\rho_{\psi}} and/or \eqn{\delta}.  Default "Percent".}  
 \item{Labels}{Logical. When \code{Labels=TRUE}, the percentage of \eqn{\rho_{\psi}} and/or \eqn{\delta} values that are equal to or larger than the midpoint value of each of the bins are displayed (on top of each bin). Default \code{FALSE}.}  
\item{PCA}{Logical. When \code{PCA=TRUE}, a plot of the PCA is provided. Default \code{TRUE}.}
\item{Good.Pretreat}{Logical. When \code{Good.Pretreat=TRUE}, a plot of \eqn{\delta} is provided. This plot is useful in the context of examinating the plausibility of finding a good pretretment predictor endpoint. For details, see Alonso et al. (submitted). Default \code{FALSE}.}
\item{EffectT0T1}{Logical. When \code{EffectT0T1=TRUE}, a plot depicting the relation between  \eqn{\rho[T0T1]} and \eqn{\rho[\psi]} is provided. Default \code{FALSE}.}
\item{Main.Good.Pretreat}{The title of the plot of \eqn{\delta}. For details, see Alonso et al. (submitted).}
  \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
  \item{col}{The color of the bins. Default \code{col <- c(8)}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{hist()}.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (submitted). Validating predictors of therapeutic success: a causal inference approach. 
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{PCA.ContCont}}
\examples{
\dontrun{
# Generate the vector of PCA.ContCont values when rho_T0S=.3, rho_T1S=.9, 
# sigma_T0T0=2, sigma_T1T1=2,sigma_SS=2, and  
# the grid of values {-1, -.99, ..., 1} is considered for the correlations
# between T0 and T1:
PCA <- PCA.ContCont(T0S=.3, T1S=.9, T0T0=2, T1T1=2, SS=2, 
T0T1=seq(-1, 1, by=.01))

# Plot the results:
plot(PCA)

# Same plot but add the percentages of PCA values that are equal to or larger 
# than the midpoint values of the bins
plot(PCA, Labels=TRUE)

# Plot of the cumulative distribution of PCA
plot(PCA, Typ="CumPerc")

}}

\keyword{Plot PCA}
\keyword{Causal-Inference framework}
\keyword{Sensitivity}
\keyword{Plausibility of a good pretreatment predictor}
\keyword{Continuous-continuous setting}