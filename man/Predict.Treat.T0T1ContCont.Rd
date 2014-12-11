\name{Predict.Treat.T0T1.ContCont}
\alias{Predict.Treat.T0T1.ContCont}
\title{
Compute the predicted treatment effect on the true endpoint of a patient based on his or her observed pretreatment predictor value in the continuous-continuous setting for a particular (single) value of \eqn{\rho_{T0T1}}.
}
\description{
This function computes the predicted \eqn{\Delta T_j} of a patient based on the pretreatment value \eqn{S_j} of a patient in the continuous-continuous setting for a particular (single) value of rho_T0T1. 
}
\usage{
Predict.Treat.T0T1.ContCont(x, S, Beta, SS, mu_S, T0T1, alpha=0.05)
}
\arguments{
  \item{x}{An object of class \code{PCA.ContCont}. See \code{\link{PCA.ContCont}}.}
  \item{S}{The observed pretreatment value \eqn{S_j} for a patient.}
  \item{Beta}{The estimated treatment effect on the true endpoint (in the validation sample).}
  \item{SS}{The estimated variance of the pretreatment predictor endpoint.}
  \item{mu_S}{The estimated mean of the surrogate endpoint (in the validation sample).}
  \item{T0T1}{The \eqn{\rho_{T0T1}} value (used to compute the variance of \eqn{\Delta T_j}|\eqn{S_j}.}
  \item{alpha}{The \eqn{\alpha}-level that is used to determine the confidence interval around \eqn{\Delta T_j}|\eqn{S_j}. Default \eqn{0.05}.}
}
\value{An object of class \code{PCA.Predict.Treat.T0T1.ContCont} with components,
\item{Pred_T}{The predicted \eqn{\Delta T_j}.}
\item{Var_Delta.T}{The variance \eqn{\sigma_{\Delta_{T}}}.}
\item{T0T1}{The correlation between the counterfactuals \eqn{T_{0}}, \eqn{T_{1}}.}
\item{CI_low}{The lower border of the \eqn{1-\alpha}\% confidence interval of \eqn{\Delta T_j}|\eqn{S_j}.}
\item{CI_high}{The upper border of the \eqn{1-\alpha}\% confidence interval of \eqn{\Delta T_j}|\eqn{S_j}.}
\item{Var_Delta.T_S}{The variance \eqn{\sigma_{\Delta_{T}}}|\eqn{S_j}.}
\item{alpha}{The \eqn{\alpha}-level that is used to determine the confidence interval of \eqn{\Delta T_j}|\eqn{S_j}.}
}

\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (submitted). Validating predictors of therapeutic success: a causal inference approach.}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{\link{PCA.ContCont}}
\examples{
# Generate the vector of PCA.ContCont values when rho_T0S=.3, rho_T1S=.9, 
# sigma_T0T0=2, sigma_T1T1=2,sigma_SS=2, and the grid of values {-1, -.99, 
# ..., 1} is considered for the correlations between T0 and T1:
PCA <- PCA.ContCont(T0S=.3, T1S=.9, T0T0=2, T1T1=2, SS=2, 
T0T1=seq(-1, 1, by=.01))

# Obtain the predicted value T for a patient who scores S = 10, using beta=5,
# SS=2, mu_S=4, assuming rho_T0T1=.6
indiv <- Predict.Treat.T0T1.ContCont(x=PCA, S=10, Beta=5, SS=2, mu_S=4, T0T1=.6)
summary(indiv)

# obtain a plot with the 95\% CI around delta T_j | S_j (assuming rho_T0T1=.6)
plot(indiv)
}
\keyword{Causal-Inference framework}
\keyword{Continuous-continuous setting}
