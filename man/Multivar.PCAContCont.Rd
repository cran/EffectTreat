\name{Multivar.PCA.ContCont}
\alias{Multivar.PCA.ContCont}
\title{Compute the multivariate predictive causal association (PCA) in the Continuous-continuous case}
\description{
The function \code{Multivar.PCA.ContCont} computes the predictive causal association (PCA) when \eqn{S} = the vector of pretreatment predictors and \eqn{T} = the True endpoint. All \eqn{S} and \eqn{T} should be continuous normally distributed endpoints. See \bold{Details} below.}
\usage{Multivar.PCA.ContCont(Sigma_TT, Sigma_TS, Sigma_SS, T0T1=seq(-1, 1, by=.01))}
\arguments{
  \item{Sigma_TT}{The variance-covariance matrix 
  \eqn{\bold{\Sigma}_{TT}=\left(\begin{array}{cc}\sigma_{T0T0} & \sigma_{T0T1} \\ \sigma_{T0T1} & \sigma_{T1T1}\end{array}\right)}.}
\item{Sigma_TS}{The matrix that contains the covariances \eqn{\sigma_{T0Sr}}, \eqn{\sigma_{T1Sr}}. For example, when there are \eqn{2} pretreatment predictors 
  \eqn{\bold{\Sigma}_{TS}=\left(\begin{array}{cc}\sigma_{T0S1} & \sigma_{T0S2} \\ \sigma_{T1S1} & \sigma_{T1S2}\end{array}\right)}.}
\item{Sigma_SS}{The variance-covariance matrix of the pretreatment predictors. For example, when there are \eqn{2} pretreatment predictors 
  \eqn{\bold{\Sigma}_{SS}=\left(\begin{array}{cc}\sigma_{S1S1} & \sigma_{S1S2} \\ \sigma_{S1S2} & \sigma_{S2S2}\end{array}\right)}.}
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals \eqn{T_0} and \eqn{T_1} that should be considered in the computation of \eqn{R^2_{\psi}}. Default \code{seq(-1, 1, by=.01)}, i.e., the values \eqn{-1}, \eqn{-0.99}, \eqn{-0.98}, \dots, \eqn{1}.}
}
\value{An object of class \code{Multivar.PCA.ContCont} with components,
\item{Total.Num.Matrices}{An object of class \code{numeric} that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.}
\item{Pos.Def}{A \code{data.frame} that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the \eqn{R^2_{\psi}} values.}
\item{PCA}{A scalar or vector that contains the PCA (\eqn{R^2_{\psi}}) value(s).}
}
\references{
Alonso, A., & Van der Elst, W. (submitted). Evaluating multivariate predictors of therapeutic success: a causal inference approach. 
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\examples{
# First specify the covariance matrices to be used 
Sigma_TT = matrix(c(177.870, NA, NA, 162.374), byrow=TRUE, nrow=2)
Sigma_TS = matrix(data = c(-45.140, -109.599, 11.290, -56.542,
-106.897, 20.490), byrow = TRUE, nrow = 2)
Sigma_SS = matrix(data=c(840.564, 73.936, -3.333, 73.936, 357.719,
-30.564, -3.333, -30.564, 95.063), byrow = TRUE, nrow = 3)

# Compute PCA
Results <- Multivar.PCA.ContCont(Sigma_TT = Sigma_TT,
Sigma_TS = Sigma_TS, Sigma_SS = Sigma_SS)

# Evaluate results
summary(Results)
plot(Results)
}

\keyword{Causal-Inference framework}
\keyword{Sensitivity}
\keyword{PCA}
\keyword{Continuous-continuous setting}
\keyword{Multivariate setting}