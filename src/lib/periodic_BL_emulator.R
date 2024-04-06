### Periodic Bayes Linear Emulator
# uses isotropic periodic covariance function with single period for all dims
# modified code from `simple_BL_emulator`

periodic_BL_emulator <- function(xD, D, xP, theta=1, sig=1, p=1, nugget=1e-3, mu=0, just_var=TRUE){
  n <- length(D); nP <- nrow(xP)
  N <- ncol(xD)
  
  Cmat <- function(X1, X2) {
    sum <- 0
    for (i in 1:N){
      dist_mat <- as.matrix(outer(X1[,i], X2[,i], FUN=function(x, y) abs(x - y)))
      sum <- sum + sin(dist_mat * pi/p)^2
    }
    return(sig^2 * exp(-4*sum/theta^2) * (1+nugget))
  }
  
  if(is.function(mu)){
    E_D <- mapply(mu, xD)
    E_B <- mapply(mu, xP)
  } else {
    E_D <- rep(mu, n)
    E_B <- rep(mu, nP)
  }
  
  Var_D <- Cmat(xD, xD)
  Cov_BD <- Cmat(xP, xD)
  Var_D_inv <- solve(Var_D)
  
  # Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x)
  Cov_BD_Var_D_inv <- Cov_BD %*% Var_D_inv # need twice, so pre-calculate
  ED_B <- E_B + Cov_BD_Var_D_inv %*% (D - E_D)
  if (just_var){
    diag_Var_B <- rep(sig^2, nP)
    
    diag_VarD_B <- diag_Var_B - apply(Cov_BD_Var_D_inv * Cov_BD, 1, sum)
    VarD_B <- diag_VarD_B  # return diagonal only
  } else {
    Var_B <- Cmat(xP, xP)
    Var_B <- Var_B + nugget*diag(nrow(Var_B)) # for stability
    
    VarD_B <- Var_B - Cov_BD_Var_D_inv %*% t(Cov_BD)
  }
  
  return(list("ExpD_f(x)"=ED_B, "VarD_f(x)"=VarD_B))
}
