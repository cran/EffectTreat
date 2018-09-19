\name{Min.R2.delta}
\alias{Min.R2.delta}
\title{Compute minimum \eqn{R^2_{\delta}} for desired prediction accuracy}
\description{
Computes the minimum \eqn{R^2_{\delta}} needed to achieve the desired prediction accuracy for the set of pretreatment predictors.}
\usage{Min.R2.delta(delta, Sigma_T0T0, Sigma_T1T1)}
\arguments{
  \item{delta}{The vector of \eqn{\delta} values to be considered.}
  \item{Sigma_T0T0}{The variance of \eqn{T} in the control treatment group.}
  \item{Sigma_T1T1}{The variance of \eqn{T} in the experimental treatment group.}
}
\references{
Alonso, A., Van der Elst, W., Luaces, P., Sanchez, L., & Molenberghs, G. (submitted). Evaluating multivariate predictors of therapeutic success: a causal inference approach.  
}
\author{
Wim Van der Elst, Ariel Alonso & Geert Molenberghs
}
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