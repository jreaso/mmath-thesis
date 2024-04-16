source("src/lib/TENSE.R")

### Simple NS Emulator V2
# no nugget term yet
simple_NS_emulator2 <- function(VxD, D, VxP, vdxD, vdxP, k_S=k_S, mu=0, theta=1, lambda2=1, sigma=1, just_var=T){
  n <- length(D); nP <- nrow(VxP)
  
  # Compute covariance matrices
  Var_D <- k_NS(VX=VxD, vdX=vdxD, k_S=k_S, theta=theta, lambda2=lambda2, sigma=sigma) # Var[xD]
  Cov_DB <- k_NS2(VX1=VxD, VX2=VxP, vdX1=vdxD, vdX2=vdxP, k_S=k_S, theta=theta, lambda2=lambda2, sigma=sigma) # Cov[xD, xP]
  
  Var_D_inv <- solve(Var_D) #chol2inv(chol(Var_D))
  Cov_BD_Var_D_inv <- t(Cov_DB) %*% Var_D_inv
  E_D <- rep(mu, n)
  E_B <- rep(mu, nP)
  
  # Perform BL expectation update
  ED_B <- E_B + Cov_BD_Var_D_inv %*% (D - E_D)
  
  if (just_var) { # compute diagonal only
    VarD_B <- rep(sigma^2, nP) - diag(Cov_BD_Var_D_inv %*% Cov_DB)
  } else {
    Var_B <- k_NS(VX=VxP, vdX=vdxP, k_S=k_S, theta=theta, lambda2=lambda2, sigma=sigma) # Var[xP]
    VarD_B <- Var_B - Cov_BD_Var_D_inv %*% t(Cov_BD)
  }
  
  return(list("ED_fP"=ED_B, "VarD_fP"=VarD_B))
}

