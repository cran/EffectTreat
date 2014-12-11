\name{plot.Predict.Treat.T0T1.ContCont}
\alias{plot Predict.Treat.T0T1.ContCont}
\alias{plot.Predict.Treat.T0T1.ContCont}
\title{
Plots the distribution of the individual causal effect based on \eqn{S} for a specific assumed correlation between the counterfactuals.}
\description{
Plots the distribution of \eqn{\Delta T_j}|\eqn{S_j} and the \eqn{1-\alpha}\% CIs for a user-requested \eqn{\rho_{T0T1}} value). The function is similar to \code{plot.Predict.Treat.ContCont}, but it is applied to an object of class \code{Predict.Treat.T0T1.ContCont} (rather than to an object of class \code{Predict.Treat.ContCont}). This object contains only one \eqn{\rho_{T0T1}} value (rather than a vector of \eqn{\rho_{T0T1}} values), and thus the plot automatically uses the considered \eqn{\rho_{T0T1}} value in the object \code{x} to compute the \eqn{1-\alpha}\% CI for \eqn{\Delta T_j}|\eqn{S_j}.   
}
\usage{
\method{plot}{Predict.Treat.T0T1.ContCont}(x, Xlab, Main, alpha=0.05, Cex.Legend=1, \dots)
}
\arguments{
  \item{x}{An object of class \code{Predict.Treat.T0T1.ContCont}. See \code{\link{Predict.Treat.T0T1.ContCont}}.}
  \item{Xlab}{The legend of the X-axis of the plot. Default "\eqn{\Delta T_j}|\eqn{S_j}".}
  \item{Main}{The title of the PCA plot. Default " ".}  
  \item{alpha}{The \eqn{\alpha} level to be used in the computation of the CIs. Default \eqn{0.05}.}
  \item{Cex.Legend}{The size of the legend of the plot. Default \eqn{1}.}
  \item{\dots}{Other arguments to be passed to the \eqn{plot()} function.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (submitted). Validating predictors of therapeutic success: a causal inference approach. 
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{Predict.Treat.T0T1.ContCont}}
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
plot(indiv, xlim=c(5, 12))
}
\keyword{Causal-Inference framework}
\keyword{Continuous-continuous setting}
\keyword{Plot}
