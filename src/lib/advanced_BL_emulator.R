### Advanced Bayes Linear Emulator - Uses Regression Surface
# uses isotropic squared exponential covariance structure

# not yet implemented just calculating diagonal of variance matrix

library(pdist)


advanced_BL_emulator <- function(xD, D, xP, mu_beta, g, Sigma_beta, theta, sig, nugget=1e-3){
  n <- length(D); nP <- nrow(xP)
  xD <- t(t(xD)/theta); xP <- t(t(xP)/theta)
  
  # Define covariance structure
  Cmat <- function(dist_mat) sig^2 * exp(-(dist_mat)^2) * (1 + nugget)
  # nugget is added to all points not just same points (when cov 0)
  
  # Define function to apply basis functions to vector
  g_x <- function(x_mat) t(apply(x_mat, 1, function(row) sapply(g, function(f) f(row))))
  
  GD <- g_x(xD) # X
  GP <- g_x(xP) # g_xP
  
  OmegaD <- Cmat(as.matrix(dist(xD)))
  CovU_D_P <- Cmat(as.matrix(pdist(xD, xP)))
  
  OmegaD_inv <- solve(OmegaD) #chol2inv(chol(OmegaD))
  Sigma_beta_inv <- solve(Sigma_beta)
  
  Gtilde <- GP - (t(CovU_D_P) %*% OmegaD_inv %*% GD) 
  
  # BL Adjustment
  VarD_beta_inv <- (t(GD) %*% OmegaD_inv %*% GD + Sigma_beta_inv)
  VarD_beta <- solve(VarD_beta_inv)
  
  ED_beta <- VarD_beta %*% (t(GD) %*% OmegaD_inv %*% D + Sigma_beta_inv %*% mu_beta)
  
  ED_U <- t(CovU_D_P) %*% OmegaD_inv %*% (D - GD %*% ED_beta)
  
  ED_B <- GP %*% ED_beta + ED_U
  
  OmegaP <- Cmat(as.matrix(dist(xP)))
  VarD_B <- OmegaP + Gtilde %*% VarD_beta %*% t(Gtilde) - t(CovU_D_P) %*% OmegaD_inv %*% CovU_D_P
  
  return(list("ExpD_f(x)"=ED_B, "VarD_f(x)"=VarD_B, "ExpD_beta"=ED_beta, "GP"=GP)) 
}
