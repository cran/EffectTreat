

# Sigma_TT=matrix(c(38.606,NA,NA,663.917),byrow=TRUE,nrow=2) 
# Sigma_TS=matrix(data=c(1101.341,2.831,7.684,2.211,-2.444,0.308,-15.403,1.916,-5.213, 
#                        -203.961,2924.364,17.195,244.002,-134.716,14.409,7.398,49.879,-3.506,94.022,-96.712), byrow=T,nrow=2)
# Sigma_SS=matrix(data=c(1016866.71,164.45,3195.37,2239.41,-225.22,128.83,-2000.46,172.40,4618.25,-3431.99, 
#                        164.448,0.827,8.980,-3.342,1.062,-0.296,9.027,0.708,1.001,-7.878, 3195.370,8.980,
#                        175.216,-42.224,12.499,-2.267,150.755,7.259,20.540,-14.477, 2239.413,-3.342,-42.224,
#                        190.735,1.060,-5.484,-23.150,3.914,-29.950,41.835, -225.226,1.0616,12.499,1.060,8.080,-2.067,37.333,-1.788,-6.060,31.833, 128.835,-0.295,-2.267,-5.484,-2.066,2.995,-28.831,-1.452,-0.541,-5.896, -2000.458,9.0267,150.755,-23.150,37.333,-28.831,886.335,-8.102,11.0241,157.235, 172.398,0.708,7.259,3.914,-1.788,-1.452,-8.102,25.885,-4.713,-62.529, 4618.248,1.001,20.540,-29.950,-6.060,-0.541,11.024,-4.713,90.247,19.882, -3431.986,-7.878,-14.477,41.835,31.833,-5.896,157.235,-62.529,19.882,1749.909), byrow=T,nrow=10)
# Results <- Multivar.PCA.ContCont(Sigma_TT = Sigma_TT, Sigma_TS = Sigma_TS, Sigma_SS = Sigma_SS)  
# plot(Results, R2_psi_g=TRUE, ylim=c(.9, 1))  
#   plot(Results)




Multivar.PCA.ContCont <- function(Sigma_TT, Sigma_TS, Sigma_SS, T0T1=seq(-1, 1, by=.01), M=NA) { 
  
  Results <- T0T1_results <- R2_psi_g_all <- NULL
  
  # If no M specified
  if (is.na(M)){
  for (i in 1: length(T0T1)) {   

    Sigma_TT_hier <- Sigma_TT
    Sigma_TT_hier[1,2] <- Sigma_TT_hier[2, 1] <- (sqrt(Sigma_TT[1,1])*sqrt(Sigma_TT[2,2])) * T0T1[i]
    
    Sigma_temp_a <- 
      rbind(cbind(Sigma_TT_hier, Sigma_TS), 
      cbind(t(Sigma_TS), Sigma_SS))
    Cor_c <- cov2cor(Sigma_temp_a)
    a_1 <- matrix(c(-1, 1), nrow=1)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)  
    
    if (Min.Eigen.Cor > 0) {
      
      PCA <- (a_1 %*% Sigma_TS %*% solve(Sigma_SS) %*% t(Sigma_TS) %*% t(a_1))/
                (a_1 %*% Sigma_TT_hier %*% t(a_1))
                 
    Results <- rbind(Results, PCA)
    T0T1_results <- rbind(T0T1_results, T0T1[i])
    }
  }
  
  if (is.null(Results)==FALSE){  
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  colnames(Results) <- "PCA"
  Pos.Def <- nrow(Results)
  }
  
  if (is.null(Results)==TRUE){  
    cat("No solutions found, try specifying another range of rho_T0T1 values. \n")
    rownames(Results) <- NULL
    Pos.Def <- 0
    T0T1_results <- NA
  }
  }
  
  
  # If M is not NA 
  if (is.na(M)==FALSE){
    
    m <- 0
    while (m < M){   
      
      Sigma_TT_hier <- Sigma_TT
      T0T1_hier <- sample(T0T1, size = 1)
      Sigma_TT_hier[1,2] <- Sigma_TT_hier[2, 1] <- (sqrt(Sigma_TT[1,1])*sqrt(Sigma_TT[2,2])) * T0T1_hier
      
      Sigma_temp_a <- 
        rbind(cbind(Sigma_TT_hier, Sigma_TS), 
              cbind(t(Sigma_TS), Sigma_SS))
      Cor_c <- cov2cor(Sigma_temp_a)
      a_1 <- matrix(c(-1, 1), nrow=1)
      Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)  
      
      if (Min.Eigen.Cor > 0) {
        m <- m+1
        PCA <- (a_1 %*% Sigma_TS %*% solve(Sigma_SS) %*% t(Sigma_TS) %*% t(a_1))/
          (a_1 %*% Sigma_TT_hier %*% t(a_1))
        
        Results <- rbind(Results, PCA)
        T0T1_results <- rbind(T0T1_results, T0T1_hier)
      }
    }
    
    if (is.null(Results)==FALSE){  
      Results <- data.frame(Results)
      rownames(Results) <- NULL
      colnames(Results) <- "PCA"
      Pos.Def <- nrow(Results)
    }
    
    if (is.null(Results)==TRUE){  
      cat("No solutions found, try specifying another range of rho_T0T1 values. \n")
      rownames(Results) <- NULL
      Pos.Def <- 0
      T0T1_results <- NA
    }
  }
  
  # Compute R^2_{psi g}
  for (s in 1: length(T0T1_results)){
  T0T1_here <- T0T1_results[s]
  Sigma_TT_here <- Sigma_TT
  Sigma_TT_here[1,2] <- Sigma_TT_here[2,1] <-  (sqrt(Sigma_TT_here[1,1])*sqrt(Sigma_TT_here[2,2])) * T0T1_here 
  R2_psi_g_here <- 1 - (det(Sigma_TT_here - Sigma_TS %*% solve(Sigma_SS) %*% t(Sigma_TS))/det(Sigma_TT_here))
  R2_psi_g_here <- c(T0T1_here, R2_psi_g_here)
  R2_psi_g_all <- rbind(R2_psi_g_all, R2_psi_g_here)
  }
  R2_psi_g <- data.frame(R2_psi_g_all, row.names = NULL)
  names(R2_psi_g) <- c("r_T0S0", "R2_psi_g")
  
  fit <- 
    list(Total.Num.Matrices=length(T0T1), Pos.Def=Pos.Def, T0T1=as.numeric(t(T0T1_results)), 
         PCA=Results$PCA, R2_psi_g=R2_psi_g, Call=match.call())
  
  class(fit) <- "Multivar.PCA.ContCont"
  fit
}
