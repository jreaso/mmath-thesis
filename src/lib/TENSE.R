###====================================
### TENSE Arbitrary Dimension
###====================================

# Stationary CF (Squared Exponential)
k_S <- function(d) exp(-d^2)

# x-dependant covariance matrix governing general Mahalanobis distance to reverse warping
Sigma_V <- function(vdX, theta=1, lambda2=1) {
  N <- length(vdX[[1]])
  theta2 <- theta^2
  
  Sigma_V_func <- function(vd) {
    r2 <- sum(vd^2)
    if (r2 == 0){ # zero derivative so no scaling
      Sigma_V <- theta2*diag(N+1)
      Sigma_V[N+1, N+1] <- lambda2
    } else {
      w_grad <- c(vd, r2)
      w_norm <- c(-vd, 1)
      Sigma_V <- ((theta2/r2) * outer(w_grad, w_grad)) + (theta2 * rbind(cbind(diag(N) - outer(vd, vd)/r2, 0), 0)) + ((lambda2/(1+r2)) * outer(w_norm, w_norm))
    }
    return(Sigma_V)
  }
  
  return(lapply(vdX, Sigma_V_func))
}


# Quadratic form for NS CF
Q <- function(x, y, Mx, My) {
  (x-y) %*% solve((Mx + My)/2) %*% (x-y)
}


# Non-Stationary TENSE Covariance Function
k_NS <- function(VX, vdX, k_S=k_S, theta=1, lambda2=1, sigma=1) {
  n <- nrow(VX)
  N <- ncol(VX) - 1
  
  Sigma_V_list <- Sigma_V(vdX, theta=theta, lambda2=lambda2)
  
  k_NS_func <- function(i, j) {
    if (i==j){
      return(1)
    } else if (j < i) {
      return(NA) # do not compute lower triangular entries
    } else {
      xi <- VX[i,]; xj <- VX[j,]
      Sigma_Vi <- Sigma_V_list[[i]]; Sigma_Vj <- Sigma_V_list[[j]]
      Qij <- Q(xi, xj, Sigma_Vi, Sigma_Vj)
      return(2^(N/2) * det(Sigma_Vi)^(1/4) * det(Sigma_Vj)^(1/4) / det(Sigma_Vi + Sigma_Vj)^(1/2) * k_S(d=sqrt(Qij)))
    }
  }
  
  Cmat <- outer(1:n, 1:n, Vectorize(k_NS_func))
  Cmat[lower.tri(Cmat)] <- t(Cmat)[lower.tri(Cmat)] # copy across diagonal
  
  return(sigma^2 * Cmat)
}
