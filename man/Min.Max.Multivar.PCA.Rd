\name{Min.Max.Multivar.PCA}
\alias{Min.Max.Multivar.PCA}
\title{Minimum and maximum values for the multivariate predictive causal association (PCA) in the continuous-continuous case}
\description{
The function \code{Min.Max.Multivar.PCA} computes the minimum and maximum values for the multivariate predictive causal association (PCA) in the continuous-continuous case.}
\usage{Min.Max.Multivar.PCA(gamma, Sigma_SS, Sigma_T0T0, Sigma_T1T1)}
\arguments{
  \item{gamma}{The vector of regression coefficients for the \eqn{S} by treatment interactions.}
  \item{Sigma_SS}{The variance-covariance matrix of the pretreatment predictors. For example, when there are \eqn{2} pretreatment predictors 
  \eqn{\bold{\Sigma}_{SS}=\left(\begin{array}{cc}\sigma_{S1S1} & \sigma_{S1S2} \\ \sigma_{S1S2} & \sigma_{S2S2}\end{array}\right)}.}
\item{Sigma_T0T0}{The variance of \eqn{T} in the control treatment group.}
\item{Sigma_T1T1}{The variance of \eqn{T} in the experimental treatment group.}
}
\references{
Alonso, A., & Van der Elst, W. (submitted). Evaluating multivariate predictors of therapeutic success: a causal inference approach. 
}
\author{
Wim Van der Elst & Ariel Alonso
}
\examples{
# Specify vector of S by treatment interaction coefficients 
gamma <- matrix(data = c(-0.006, -0.002, 0.045), ncol=1)
# Specify variances
Sigma_SS = matrix(data=c(882.352, 49.234, 6.420,
49.234, 411.964, -26.205, 6.420, -26.205, 95.400),
byrow = TRUE, nrow = 3)
Sigma_T0T0 <- 82.274
Sigma_T1T1 <- 96.386

# Compute min and max PCA
Min.Max.Multivar.PCA(gamma=gamma, Sigma_SS=Sigma_SS, 
Sigma_T0T0=Sigma_T0T0, Sigma_T1T1=Sigma_T1T1)
}
\keyword{min, max PCA}
\keyword{Continuous-continuous setting}
\keyword{Multivariate setting}