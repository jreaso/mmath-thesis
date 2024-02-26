library('SimDesign')  # For sampling from MV Normal

# =============================================
# === 1D GP WITH SQUARED EXPONENTIAL KERNEL
# =============================================

# Squared Exponential Kernel
kernel <- function(dist_matrix) exp(-(dist_matrix)^2)

# 1D GP
simple_1d_gp <- function(xP, xD, D, sigma=1, theta=1){
  n <- length(D); nP <- length(xP)
  
  xP <- t(t(xP)/theta); xD <- t(t(xD)/theta)
  
  E_D <- rep(0, n); E_fx <- rep(0, nP)
  
  Cov_fx_fxdash <- function(dist_matrix) {
    sigma^2 * kernel(dist_matrix)
  }
  
  Var_D <- Cov_fx_fxdash(as.matrix(dist(xD)))
  Var_D_inv <- chol2inv(chol(Var_D))
  
  Var_fx <- Cov_fx_fxdash(as.matrix(dist(xP)))
  
  Cov_fx_D <- Cov_fx_fxdash(as.matrix(abs(sapply(xD, function(y) (outer(xP, y, "-"))))))
  
  # Update
  ED_fx <- E_fx + Cov_fx_D %*% Var_D_inv %*% (D - E_D)
  VarD_fx <- Var_fx - Cov_fx_D %*% Var_D_inv %*% t(Cov_fx_D)
  
  return(list("Exp"=ED_fx, "Var"=VarD_fx))
}


# =============================================
# === GP PLOTTING FUNCTION
# =============================================

plot_gp <- function(gp_exp, gp_var, xP, xD=NULL, D=NULL, samples=NULL, fxP=NULL, ylim=NULL){
  int_col <- rgb(.9,.9,.9)
  sample_col <- rgb(.7,.7,.7)
  
  # Initialise the plot window
  plot(1, 1, xlim=c(0,1), ylim=ylim, type="n", xlab="", ylab="")
  
  # Plot the prediction interval
  upper <- gp_exp + 2*sqrt(gp_var); lower <- gp_exp - 2*sqrt(gp_var)
  polygon(c(xP, rev(xP)), c(upper, rev(lower)), col=int_col, border=NA)
  
  # Plot samples
  if (!is.null(samples)) {
    for (i in 1:dim(samples)[1]){
      lines(xP, as.vector(samples[i,]), col=sample_col, lwd=1)
    }
  }
  
  # Plot true function
  if (!is.null(fxP)) {
    lines(xP, fxP, col='red', lwd=2.5)
  }
  
  # Plot the gp expectation
  lines(xP, gp_exp, col=1, lwd=2, lty=2.5)
  
  # Plot the evaluation points
  if (!is.null(xD)){
    points(xD, D, pch=21, col=1, bg=1, cex=1)
  }
  
  box()
}


# =============================================
# === PRIOR PLOT
# =============================================

f <- function(x) sin(2*pi*x+5) + cos(3*pi*x+1)

xP <- seq(-.05, 1.05, len=1001)

sigma <- 1
theta <- .2

gp_exp_prior <- rep(0, length(xP))
gp_var_prior_mat <- sigma^2 * kernel(as.matrix(dist( t(t(xP)/theta) )))  # 1D so can scale theta after dist
gp_var_prior <- diag(gp_var_prior_mat)

# Sample from posterior
set.seed(42)
n_samples <- 12
prior_samples <- rmvnorm(n=n_samples, mean=gp_exp_prior, sigma=gp_var_prior_mat)

pdf(file="figures/1d-gp.pdf", width=6, height=5)
par(mar = c(4, 4, 2, 2))
plot_gp(gp_exp=gp_exp_prior, gp_var=gp_var_prior, xP=xP, 
        samples=prior_samples, fxP=f(xP), ylim=c(-3.5, 3.5))
title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)


# =============================================
# === POSTERIOR GP ON 4 POINTS
# =============================================

xD <- seq(.125,.875, .25)
D <- f(xD)

gp_out <- simple_1d_gp(xP, xD, D, sigma=sigma, theta=theta)

gp_exp <- as.vector(gp_out[["Exp"]])
gp_var_mat <- as.matrix(gp_out[["Var"]])
gp_var <- diag(gp_var_mat)

# Sample from posterior
set.seed(42)
samples <- rmvnorm(n=n_samples, mean=gp_exp, sigma=gp_var_mat)

plot_gp(gp_exp=gp_exp, gp_var=gp_var, xP=xP, xD=xD, D=D, samples=samples, fxP=f(xP), ylim=c(-3.5, 3.5))
title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)
dev.off()

