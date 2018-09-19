\name{plot Min.R2.delta}
\alias{plot Min.R2.delta}
\alias{plot.Min.R2.delta}
\title{
Plot \eqn{R^2_{\delta}} as a function of \eqn{\delta}.
}
\description{
This function plots \eqn{R^2_{\delta}} as a function of \eqn{\delta} (in the multivariate case).}
\usage{
\method{plot}{Min.R2.delta}(x, Ylab, Main="", Ylim=c(0, 1), \dots)
}
\arguments{
  \item{x}{An object of class \code{plot.Min.R2.delta}. See \code{\link{Min.R2.delta}}.}
\item{Ylab}{The legend of the Y-axis of the PCA plot. Default \eqn{R^2_{\delta}}.}
 \item{Main}{The title of the plot. Default " " (no title).}  
 \item{Ylim}{The limits of the Y-axis.  Default \code{Ylim=c(0, 1)}.}  
  \item{\dots}{Extra graphical parameters to be passed to \code{plot()}.}
}
\references{
Alonso, A., Van der Elst, W., Luaces, P., Sanchez, L., & Molenberghs, G. (submitted). Evaluating multivariate predictors of therapeutic success: a causal inference approach.  
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{Min.R2.delta}}
\examples{
Fit <- Min.R2.delta(delta = seq(from = 0, to = 250, by=50), 
  Sigma_T0T0 = 38.606, Sigma_T1T1 = 663.917)

# Explore the results
summary(Fit)
plot(Fit)
}
\keyword{Prediction accuracy}
\keyword{Continuous-continuous setting}
\keyword{Multivariate setting}
\keyword{Multivariate predictors}
