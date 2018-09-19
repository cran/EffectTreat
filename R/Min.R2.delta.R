#Fit <- Min.R2.delta(delta = seq(from = 0, to = 250, by=50), Sigma_T0T0 = 38.606, Sigma_T1T1 = 663.917)
#summary(Fit)
#plot(Fit)

Min.R2.delta <- function(delta, Sigma_T0T0, Sigma_T1T1){ 

R2_delta <- 1 - ((delta / ((sqrt(Sigma_T0T0) + sqrt(Sigma_T1T1))**2)))
R2_delta[R2_delta<0] <- 0
    
fit <- 
    list(R2_delta=R2_delta, delta=delta, Call=match.call())
  
  class(fit) <- "Min.R2.delta"
  fit
}




summary.Min.R2.delta <- function(object, ..., Object){
  if (missing(Object)){Object <- object} 
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n# R2_delta for different delta's:")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  temp <- data.frame(Object$delta, Object$R2_delta)
  names(temp) <- c("delta", "R^2_delta")
  print(temp)
}


plot.Min.R2.delta <- function(x, Ylab, Main="", Ylim=c(0, 1), ...){ 
  
  Object <- x 
  
    if (missing(Ylab)) {Ylab <- expression(paste(R[delta]^2))}
    par(mar=c(4,4.5,4,2) + 0.1) 
    plot(x=Object$delta, y=Object$R2_delta, type="l", xlab=expression(delta), ylab=Ylab, main=Main, lwd=2, col=1,
         ylim=Ylim, ...)  
    par(mar=c(5,4,4,2) + 0.1) 
    
}
