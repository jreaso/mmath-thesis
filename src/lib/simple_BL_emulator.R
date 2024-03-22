### Simple Bayes Linear Emulator
# uses isotropic squared exponential covariance structure

library(pdist)

simple_BL_emulator <- function(xD, D, xP, theta=1, sig=1, nugget=1e-3, mu=0, just_var=TRUE){
  n <- length(D); nP <- nrow(xP)
  
  Cmat <- function(dist_mat) sig^2 * exp(-(dist_mat)^2) * (1+nugget)
  
  if(is.function(mu)){
    E_D <- mapply(mu, xD)
    E_B <- mapply(mu, xP)
  } else {
    E_D <- rep(mu, n)
    E_B <- rep(mu, nP)
  }
  
  xD <- t(t(xD)/theta); xP <- t(t(xP)/theta)
  
  Var_D <- Cmat(as.matrix(dist(xD)))
  Var_D <- Var_D + nugget*diag(nrow(Var_D)) # for stability
  Cov_BD <- Cmat(as.matrix(pdist(xP, xD)))
  Var_D_inv <- chol2inv(chol(Var_D))
  
  # Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x)
  Cov_BD_Var_D_inv <- Cov_BD %*% Var_D_inv # need twice, so precalculate
  ED_B <- E_B + Cov_BD_Var_D_inv %*% (D - E_D)
  if (just_var){
    diag_Var_B <- rep(sig^2, nP)
    
    diag_VarD_B <- diag_Var_B - apply(Cov_BD_Var_D_inv * Cov_BD, 1, sum)
    VarD_B <- diag_VarD_B  # return diagonal only
  } else {
    Var_B <- Cmat(as.matrix(dist(xP)))
    Var_B <- Var_B + nugget*diag(nrow(Var_B)) # for stability
    
    VarD_B <- Var_B - Cov_BD_Var_D_inv %*% t(Cov_BD)
  }
  
  return(list("ExpD_f(x)"=ED_B, "VarD_f(x)"=VarD_B))
}
