source("src/lib/TENSE.R")
  
### Advanced NS Emulator
# no nugget term yet and only computes variances, not variance matrix
advanced_NS_emulator <- function(VxD, D, VxP, vdxD, vdxP, k_S=k_S, mu_beta, g, Sigma_beta, theta=1, lambda2=1, sigma=1){ #nugget=1e-3
  n <- length(D); nP <- nrow(xP)
  xD <- t(t(xD)/theta); xP <- t(t(xP)/theta)
  
  # Compute covariance matrices
  OmegaD <- k_NS(VX=VxD, vdX=vdxD, k_S=k_S, theta=theta, lambda2=lambda2, sigma=sigma) # Var[xD]
  OmegaP <- k_NS(VX=VxP, vdX=vdxP, k_S=k_S, theta=theta, lambda2=lambda2, sigma=sigma) # Var[xP]
  CovU_D_P <- k_NS2(VX1=VxD, VX2=VxP, vdX1=vdxD, vdX2=vdxP, k_S=k_S, theta=theta, lambda2=lambda2, sigma=sigma) # Cov[xD, xP]
  
  # Apply basis functions to points
  g_x <- function(x_mat) t(apply(x_mat, 1, function(row) sapply(g, function(f) f(row))))
  GD <- g_x(xD)
  GP <- g_x(xP)
  
  OmegaD_inv <- solve(OmegaD) #chol2inv(chol(OmegaD))
  Sigma_beta_inv <- solve(Sigma_beta)
  
  Gtilde <- GP - (t(CovU_D_P) %*% OmegaD_inv %*% GD) 
  
  # BL Adjustment
  VarD_beta_inv <- (t(GD) %*% OmegaD_inv %*% GD + Sigma_beta_inv)
  VarD_beta <- solve(VarD_beta_inv)
  
  ED_beta <- VarD_beta %*% (t(GD) %*% OmegaD_inv %*% D + Sigma_beta_inv %*% mu_beta)
  ED_U <- t(CovU_D_P) %*% OmegaD_inv %*% (D - GD %*% ED_beta)
  ED_B <- GP %*% ED_beta + ED_U
  
  VarD_B <- OmegaP + Gtilde %*% VarD_beta %*% t(Gtilde) - t(CovU_D_P) %*% OmegaD_inv %*% CovU_D_P
  
  return(list("ED_fP"=ED_B, "VarD_fP"=VarD_B, "ED_beta"=ED_beta, "GP"=GP)) 
}
  