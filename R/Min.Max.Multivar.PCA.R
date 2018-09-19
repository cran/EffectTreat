Min.Max.Multivar.PCA <- function(gamma, Sigma_SS, Sigma_T0T0, Sigma_T1T1){ 

  Max.PCA <- min(t(gamma) %*% Sigma_SS %*% gamma / ( (sqrt(Sigma_T0T0) - sqrt(Sigma_T1T1)   )**2), 1) 
  Min.PCA <- max(t(gamma) %*% Sigma_SS %*% gamma / ( (sqrt(Sigma_T0T0) + sqrt(Sigma_T1T1)   )**2), 0)
  
  cat("\n \nMin PCA: ", Min.PCA)
  cat("\n \nMax PCA: ", Max.PCA, "\n \n")  
  
fit <- 
    list(Call=match.call())
  
  class(fit) <- "Min.Max.Multivar.PCA"
  fit
}
