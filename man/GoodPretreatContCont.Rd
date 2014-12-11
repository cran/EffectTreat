\name{GoodPretreatContCont}
\alias{GoodPretreatContCont}
\title{Examine the plausibility of finding a good pretreatment predictor in the Continuous-continuous case}

\description{
The function \code{GoodPretreatContCont} examines the plausibility of finding a good pretreatment predictor in the continuous-continuous setting. For details, see Alonso et al. (submitted).}

\usage{GoodPretreatContCont(T0T0, T1T1, Delta, T0T1=seq(from=0, to=1, by=.01))}
\arguments{
  \item{T0T0}{A scalar that specifies the variance of the true endpoint in the control treatment condition.}
  \item{T1T1}{A scalar that specifies the variance of the true endpoint in the experimental treatment condition.}
  \item{Delta}{A scalar that specifies an upper bound for the prediction mean squared error when predicting the individual causal effect of the treatment on the true endpoint based on the pretreatment predictor.}
  \item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals \eqn{T_0} and \eqn{T_1} that should be considered in the computation of \eqn{\rho_{min}^{2}}. Default \code{seq(0, 1, by=.01)}, i.e., the values \eqn{0}, \eqn{0.01}, \eqn{0.02}, \dots, \eqn{1}.}
  }


\value{An object of class \code{GoodPretreatContCont} with components, 
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1 that were considered (i.e., \eqn{\rho(_{T_{0},T_{1}})}).}
\item{Sigma.Delta.T}{A scalar or vector that contains the standard deviations of the individual causal treatment effects on the true endpoint as a function of \eqn{\rho(_{T_{0},T_{1})}}.}
\item{Rho2.Min}{A scalar or vector that contains the \eqn{\rho_{min}^{2}} values as a function of \eqn{\rho(_{T_{0},T_{1}})}.}
}

\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (submitted). Validating predictors of therapeutic success: a causal inference approach. 
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{
\code{\link{PCA.ContCont}}
}

\examples{
# Assess the plausibility of finding a good pretreatment predictor when
# sigma_T0T0 = sigma_T1T1 = 8 and Delta = 1
MinPred <- GoodPretreatContCont(T0T0 = 8, T1T1 = 8, Delta = 1)
summary(MinPred)
\donttest{plot(MinPred)}
}
\keyword{Plausibility of a good pretreatment predictor}
\keyword{Continuous-continuous setting}