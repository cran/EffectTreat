\name{Predict.Treat.Multivar.ContCont}
\alias{Predict.Treat.Multivar.ContCont}
\title{
Compute the predicted treatment effect on the true endpoint of a patient based on his or her observed vector of pretreatment predictor values in the continuous-continuous setting
}
\description{
This function computes the predicted \eqn{\Delta T_j} of a patient based on the vector of pretreatment values \eqn{\bold{S}_j} of a patient in the continuous-continuous setting. 
}
\usage{
Predict.Treat.Multivar.ContCont(Sigma_TT, Sigma_TS, Sigma_SS, Beta, 
S, mu_S, T0T1=seq(-1, 1, by=.01))
}
\arguments{
  \item{Sigma_TT}{The variance-covariance matrix 
  \eqn{\bold{\Sigma}_{TT}=\left(\begin{array}{cc}\sigma_{T0T0} & \sigma_{T0T1} \\ \sigma_{T0T1} & \sigma_{T1T1}\end{array}\right)}.}
\item{Sigma_TS}{The matrix that contains the covariances \eqn{\sigma_{T0Sr}}, \eqn{\sigma_{T1Sr}}. For example, when there are \eqn{2} pretreatment predictors 
  \eqn{\bold{\Sigma}_{TS}=\left(\begin{array}{cc}\sigma_{T0S1} & \sigma_{T0S2} \\ \sigma_{T1S1} & \sigma_{T1S2}\end{array}\right)}.}
\item{Sigma_SS}{The variance-covariance matrix of the pretreatment predictors. For example, when there are \eqn{2} pretreatment predictors 
  \eqn{\bold{\Sigma}_{SS}=\left(\begin{array}{cc}\sigma_{S1S1} & \sigma_{S1S2} \\ \sigma_{S1S2} & \sigma_{S2S2}\end{array}\right)}.}
  \item{Beta}{The estimated treatment effect on the true endpoint (in the validation sample).}
  \item{S}{The vector of observed pretreatment values \eqn{\bold{S}_j} for a patient.}
  \item{mu_S}{The vector of estimated means of the pretreatment predictor (in the validation sample).}
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals \eqn{T_0} and \eqn{T_1} that should be considered in the computation of \eqn{\rho_{\psi}}. Default \code{seq(-1, 1, by=.01)}, i.e., the values \eqn{-1}, \eqn{-0.99}, \eqn{-0.98}, \dots, \eqn{1}.}
}
\value{An object of class \code{PCA.Predict.Treat.Multivar.ContCont} with components,
\item{Pred_T}{The predicted \eqn{\Delta T_j}.}
\item{Var_Delta.T_S}{The variance \eqn{\sigma_{\Delta_{T}}}|\eqn{S_j}.}
\item{T0T1}{The correlation between the counterfactuals \eqn{T_{0}}, \eqn{T_{1}}.}
}
\references{
Alonso, A., & Van der Elst, W. (submitted). Evaluating multivariate predictors of therapeutic success: a causal inference approach. }

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{\link{PCA.ContCont}, \link{Multivar.PCA.ContCont}}
\examples{
# Specify the covariance matrices to be used 
Sigma_TT = matrix(c(177.870, NA, NA, 162.374), byrow=TRUE, nrow=2)
Sigma_TS = matrix(data = c(-45.140, -109.599, 11.290, -56.542,
-106.897, 20.490), byrow = TRUE, nrow = 2)
Sigma_SS = matrix(data=c(840.564, 73.936, -3.333, 73.936, 357.719,
-30.564, -3.333, -30.564, 95.063), byrow = TRUE, nrow = 3)

# Specify treatment effect (Beta), means of vector S (mu_s), and 
# observed pretreatment variable values for patient (S)
Beta <- -0.9581 # treatment effect
mu_S = matrix(c(66.8149, 84.8393, 25.1939), nrow=3) #means S_1--S_3
S = matrix(c(90, 180, 30), nrow=3) # S_1--S_3 values for a patient

# predict Delta_T based on S
Pred_S <- Predict.Treat.Multivar.ContCont(Sigma_TT=Sigma_TT, Sigma_TS=Sigma_TS,
Sigma_SS=Sigma_SS, Beta=Beta, S=S, mu_S=mu_S, T0T1=seq(-1, 1, by=.01))

# Explore results
summary(Pred_S)
plot(Pred_S)
}
\keyword{Causal-Inference framework}
\keyword{Continuous-continuous setting}
\keyword{Multivariate setting}