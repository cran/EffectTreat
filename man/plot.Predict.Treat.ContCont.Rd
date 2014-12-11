\name{plot.Predict.Treat.ContCont}
\alias{plot Predict.Treat.ContCont}
\alias{plot.Predict.Treat.ContCont}
\title{
Plots the distribution of the individual causal effect based on \eqn{S}. 
}
\description{
Plots the distribution of \eqn{\Delta T_j}|\eqn{S_j} and the \eqn{1-\alpha}\% CIs for the mean and median \eqn{\rho_{T0T1}} values (and optionally, for other user-requested \eqn{\rho_{T0T1}} values).
}
\usage{
\method{plot}{Predict.Treat.ContCont}(x, Xlab, Main, Mean.T0T1=FALSE, Median.T0T1=TRUE,  
Specific.T0T1="none", alpha=0.05, Cex.Legend=1, \dots)
}
\arguments{
  \item{x}{An object of class \code{Predict.Treat.ContCont}. See \code{\link{Predict.Treat.ContCont}}.}
  \item{Xlab}{The legend of the X-axis of the plot. Default "\eqn{\Delta T_j}|\eqn{S_j}".}
  \item{Main}{The title of the PCA plot. Default " ".}  
  \item{Mean.T0T1}{Logical. When \code{Mean.T0T1=TRUE}, the \eqn{1-\alpha}\% CI for the mean \eqn{\rho_{T0T1}} value (i.e., the mean of all valid \eqn{\rho_{T0T1}} values in \code{x}) is shown. Default \code{FALSE}.}  
  \item{Median.T0T1}{Logical. When \code{Median.T0T1=TRUE}, the \eqn{1-\alpha}\% CI for the median \eqn{\rho_{T0T1}} value is shown. Default \code{TRUE}.} 
  \item{Specific.T0T1}{Optional. A scalar that specifies a particular value \eqn{\rho_{T0T1}} for which the \eqn{1-\alpha}\% CI is shown. Default \code{"none"}.}
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
\seealso{\link{Predict.Treat.ContCont}}
\examples{
# Generate the vector of PCA.ContCont values when rho_T0S=.3, rho_T1S=.9, 
# sigma_T0T0=2, sigma_T1T1=2,sigma_SS=2, and the grid of values {-1, -.99, 
# ..., 1} is considered for the correlations between T0 and T1:
PCA <- PCA.ContCont(T0S=.3, T1S=.9, T0T0=2, T1T1=2, SS=2, 
T0T1=seq(-1, 1, by=.01))

# Obtain the predicted value T for a patient who scores S = 10, using beta=5,
# SS=2, mu_S=4
Predict <- Predict.Treat.ContCont(x=PCA, S=10, Beta=5, SS=2, mu_S=4)

# examine the results
summary(Predict)

# plot Delta_T_j given S_T and 95\% CI  based on 
# the mean value of the valid rho_T0T1 results  
plot(Predict, Mean.T0T1=TRUE, Median.T0T1=FALSE, 
xlim=c(4, 13))

# plot Delta_T_j given S_T and 99\% CI  using 
#  rho_T0T1=.8   
plot(Predict, Mean.T0T1=FALSE, Median.T0T1=FALSE, 
Specific.T0T1=.6, alpha=0.01, xlim=c(4, 13))
}
\keyword{Causal-Inference framework}
\keyword{Continuous-continuous setting}
\keyword{Plot}
