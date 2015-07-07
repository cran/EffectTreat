\name{CausalPCA.ContCont}
\alias{CausalPCA.ContCont}
\title{
Show a causal diagram of the median correlation between the counterfactuals in the continuous-continuous setting
}
\description{
This function provides a diagram that depicts the estimable correlations \eqn{\rho(_{T_0, S})} and \eqn{\rho(_{T_1, S})}, and median of the correlation \eqn{\rho(_{T_0, T_1})} for a specified range of values of the predictive causal association (PCA; \eqn{\rho_{\psi}}). 
}
\usage{
CausalPCA.ContCont(x, Min=-1, Max=1, Cex.Letters=3, Cex.Corrs=2, 
Lines.Rel.Width=TRUE, Col.Pos.Neg=TRUE)
}
\arguments{
  \item{x}{An object of class \code{PCA.ContCont}. See \code{\link{PCA.ContCont}}.}
  \item{Min}{The minimum values of the PCA that should be considered. Default=\eqn{-1}.}
  \item{Max}{The maximum values of the PCA that should be considered. Default=\eqn{1}.}
  \item{Cex.Letters}{The size of the symbols for \eqn{S}, \eqn{T_{0}}, and \eqn{T_{1}} in the diagram. Default=\eqn{3}.}
  \item{Cex.Corrs}{The size of the text depicting the (median) correlations in the diagram. Default=\eqn{2}.}
  \item{Lines.Rel.Width}{Logical. When \code{Lines.Rel.Width=TRUE}, the widths of the lines that represent the correlations in the diagram are relative to the size of the correlations (i.e., a smaller line is used for correlations closer to zero whereas a thicker line is used for (absolute) correlations closer to \eqn{1}). When \code{Lines.Rel.Width=FALSE}, the width of all lines representing the correlations between the counterfactuals is identical. Default=\code{TRUE}.}
  \item{Col.Pos.Neg}{Logical. When \code{Col.Pos.Neg=TRUE}, the color of the lines that represent the correlations in the diagram is red for negative correlations and black for positive ones. When \code{Col.Pos.Neg=FALSE}, all lines are in black. Default=\code{TRUE}.}
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

# Obtain causal diagram for PCA score range [-1; 1]:
CausalPCA.ContCont(PCA,  Min=-1, Max=1)

# Obtain causal diagram for PCA score range [0.5; 1]:
CausalPCA.ContCont(PCA,  Min=0.5, Max=1)
}
\keyword{Causal-Inference framework}
\keyword{Continuous-continuous setting}
