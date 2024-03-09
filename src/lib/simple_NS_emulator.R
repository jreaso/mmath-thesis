### Bayes Linear Non-stationary Emulator
simple_NS_emulator <- function(xD, D, xP, sig=1, nugget=1e-3, mu=0, Sigma, just_var=TRUE){
  n <- nrow(xD); nP <- nrow(xP)
  
  E_B <- rep(mu, nP)
  E_D <- rep(mu, n)
  
  # VECTORISE k_NS AND REMOVE THESE FOR LOOPS FOR BETTER EFFICIENCY
  k_nug <- function(x, y) (1-nugget) * k_NS(x=x, y=y, sig=sig, Sigma=Sigma)
  
  Var_D <- diag(sig^2/2, n)  # half sig^2 as gets doubled when transpose added
  # only calc upper diagonal
  for(i in 1:(n-1)) for(j in (i+1):n)  Var_D[i,j] <- (1-nugget)*k_NS(x=xD[i, ], y=xD[j, ], sig=sig, Sigma=Sigma)
  Var_D <- Var_D + t(Var_D)  # add transpose to form full matrix, now with correct sig^2 diagonal
  
  Cov_BD <- matrix(0, nrow=nP, ncol=n)
  for(i in 1:nP) for(j in 1:n)  Cov_BD[i,j] <- (1-nugget)*k_NS(x=xP[i, ], y=xD[j, ], sig=sig, Sigma=Sigma)
  
  
  Var_D_inv <- chol2inv(chol(Var_D))
  Cov_BD_Var_D_inv <- Cov_BD %*% Var_D_inv  # need twice, so precalculate
  
  # Perform BL expectation update
  ED_B <- E_B + Cov_BD_Var_D_inv %*% (D - E_D)
  
  if (just_var) {
    diag_Var_B <- rep(sig^2, nP)
    
    diag_VarD_B <- diag_Var_B - apply(Cov_BD_Var_D_inv * Cov_BD, 1, sum)
    VarD_B <- diag_VarD_B  # return diagonal only
  } else {
    Var_B <- diag(sig^2/2, nP)
    for(i in 1:(nP-1)) for(j in (i+1):nP)  Var_D[i,j] <- (1-nugget)*k_NS(x=xP[i, ], y=xP[j, ], sig=sig, Sigma=Sigma)
    Var_B <- Var_B + t(Var_B)  # add transpose to form full matrix, now with correct sig^2 diagonal
    
    VarD_B <- Var_B - Cov_BD_Var_D_inv %*% t(Cov_BD)
  }
  
  return(list("ExpD_f(x)"=ED_B, "VarD_f(x)"=VarD_B))
}