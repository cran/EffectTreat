\name{PCA.ContCont}
\alias{PCA.ContCont}
\title{Compute the predictive causal association (PCA) in the Continuous-continuous case}
\description{
The function \code{PCA.ContCont} compute the predictive causal association (PCA) when \eqn{S}=pretreatment predictor and \eqn{T}=True endpoint are continuous normally distributed endpoints. See \bold{Details} below.}
\usage{PCA.ContCont(T0S, T1S, T0T0=1, T1T1=1, SS=1, T0T1=seq(-1, 1, by=.01))}
\arguments{
  \item{T0S}{A scalar or vector that specifies the correlation(s) between the pretreatment predictor and the true endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{\psi}}.}
  \item{T1S}{A scalar or vector that specifies the correlation(s) between the pretreatment predictor and the true endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{\psi}}.}
  \item{T0T0}{A scalar that specifies the variance of the true endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{\psi}}. Default 1.}
  \item{T1T1}{A scalar that specifies the variance of the true endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{\psi}}. Default 1.}
  \item{SS}{A scalar that specifies the variance of the pretreatment predictor endpoint. Default 1.}
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals \eqn{T_0} and \eqn{T_1} that should be considered in the computation of \eqn{\rho_{\psi}}. Default \code{seq(-1, 1, by=.01)}, i.e., the values \eqn{-1}, \eqn{-0.99}, \eqn{-0.98}, \dots, \eqn{1}.}
}
\details{
Based on the causal-inference framework, it is assumed that each subject \emph{j} has two counterfactuals (or potential outcomes), i.e., \eqn{T_{0j}} and \eqn{T_{1j}} (the counterfactuals for the true endpoint (\eqn{T}) under the control (\eqn{Z=0}) and the experimental (\eqn{Z=1}) treatments of subject \emph{j}, respectively). The individual causal effects of \eqn{Z} on \eqn{T} for a given subject \emph{j} is then defined as \eqn{\Delta_{T_{j}}=T_{1j}-T_{0j}}.

Surrogacy can be quantified as the correlation between the individual causal effect of \eqn{Z} on \eqn{T} and \eqn{S_{j}} (the pretreatment predictor) (for details, see Alonso et al., submitted):

\deqn{\rho_{\psi}=\frac{\sqrt{\sigma_{T1T1}}\rho_{T1S}-\sqrt{\sigma_{T0T0}}\rho_{T0S}}{\sqrt{\sigma_{T0T0}+\sigma_{T1T1}-2\sqrt{\sigma_{T0T0}\sigma_{T1T1}}}\rho_{T0T1}},}

where the correlation \eqn{\rho_{T_{0}T_{1}}} is not estimable. It is thus warranted to conduct a sensitivity analysis (by considering vectors of possible values for the correlations between the counterfactuals -- rather than point estimates). 

When the user specifies a vector of values that should be considered for \eqn{\rho_{T_{0}T_{1}}} in the above expression, the function \code{PCA.ContCont} constructs all possible matrices that can be formed as based on these values and the estimable quantities \eqn{\rho_{T_{0}S}}, \eqn{\rho_{T_{1}S}}, identifies the matrices that are positive definite (i.e., valid correlation matrices), and computes \eqn{\rho_{\psi}} for each of these matrices. The obtained vector of \eqn{\rho_{\psi}} values can subsequently be used to e.g., conduct a sensitivity analysis.   

\bold{Notes}

A single \eqn{\rho_{\psi}} value is obtained when all correlations in the function call are scalars.

}
\value{An object of class \code{PCA.ContCont} with components,
\item{Total.Num.Matrices}{An object of class \code{numeric} that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.}
\item{Pos.Def}{A \code{data.frame} that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the \eqn{\rho_{\psi}} values.}
\item{PCA}{A scalar or vector that contains the PCA (\eqn{\rho_{\psi}}) value(s).}
\item{GoodSurr}{A \code{data.frame} that contains the PCA (\eqn{\rho_{\psi}}), \eqn{\sigma_{\psi_{T}}}, and \eqn{\delta}.}  
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (submitted). Validating predictors of therapeutic success: a causal inference approach. 
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\examples{

# Based on the example dataset
    # load data in memory
data(Example.Data)
    # compute corr(S, T) in control treatment, gives .77
cor(Example.Data$S[Example.Data$Treat==-1], 
Example.Data$T[Example.Data$Treat==-1])
   # compute corr(S, T) in experimental treatment, gives .71
cor(Example.Data$S[Example.Data$Treat==1], 
Example.Data$T[Example.Data$Treat==1])
   # compute var T in control treatment, gives 263.99 
var(Example.Data$T[Example.Data$Treat==-1])
   # compute var T in experimental treatment, gives 230.64  
var(Example.Data$T[Example.Data$Treat==1])
   # compute var S, gives 163.65   
var(Example.Data$S)

# Generate the vector of PCA.ContCont values using these estimates 
# and the grid of values {-1, -.99, ..., 1} for the correlations
# between T0 and T1:
PCA <- PCA.ContCont(T0S=.77, T1S=.71, T0T0=263.99, T1T1=230.65, 
                    SS=163.65, T0T1=seq(-1, 1, by=.01))

# Examine and plot the vector of generated PCA values:
summary(PCA)
\donttest{plot(PCA)}

# Other example

# Generate the vector of PCA.ContCont values when rho_T0S=.3, rho_T1S=.9, 
# sigma_T0T0=2, sigma_T1T1=2,sigma_SS=2, and  
# the grid of values {-1, -.99, ..., 1} is considered for the correlations
# between T0 and T1:
PCA <- PCA.ContCont(T0S=.3, T1S=.9, T0T0=2, T1T1=2, SS=2, 
T0T1=seq(-1, 1, by=.01))

# Examine and plot the vector of generated PCA values:
summary(PCA)
\donttest{plot(PCA)}

# Obtain the positive definite matrices than can be formed as based on the 
# specified (vectors) of the correlations (these matrices are used to 
# compute the PCA values)
PCA$Pos.Def
}

\keyword{Causal-Inference framework}
\keyword{Sensitivity}
\keyword{PCA}
\keyword{Continuous-continuous setting}