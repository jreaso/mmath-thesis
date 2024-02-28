library(SimDesign)


# =============================================
# === COVARIANCE STRUCTURES
# =============================================

sigma <- 1; theta <- 0.25

# (Isotropic) Squared exponential kernel
se_kernel <- function(dist_mat) exp(-(dist_mat)^2)

# Matérn kernel
matern_kernel <- function(dist_mat, nu=1.5) {
  A = 2*sqrt(nu) * dist_mat
  ifelse(dist_mat==0, 1, 2^(1-nu) / gamma(nu) * (A)^nu * besselK(A, nu))
}

# Exponential kernel
exp_kernel <- function(dist_mat) exp(-dist_mat)


# =============================================
# === 1D PLOTS
# =============================================

n <- 1001
x_seq <- seq(-0.1, 1.1, l=n)

dmat <- as.matrix(dist(x_seq/theta))

se_cov_mat <- sigma^2 * se_kernel(dmat)
matern_cov_mat <- sigma^2 * matern_kernel(dmat)
exp_cov_mat <- sigma^2 * exp_kernel(dmat)

set.seed(42)
se_samples <- rmvnorm(n=12, mean=rep(0, n), sigma=se_cov_mat)
matern_samples <- rmvnorm(n=12, mean=rep(0, n), sigma=matern_cov_mat)
exp_samples <- rmvnorm(n=5, mean=rep(0, n), sigma=exp_cov_mat)


sample_col <- rgb(.7,.7,.7)

# 1D Sample Plotting Function
plot_samples <- function(x_seq, samples, ylim=NULL, lwd=1){
  # Initialise the plot window
  plot(1, 1, xlim=c(0,1), ylim=ylim, type="n", xlab="", ylab="")
  
  # Plot samples
  for (i in 1:dim(samples)[1]){
    lines(x_seq, as.vector(samples[i,]), col=sample_col, lwd=lwd)
  }
}


pdf(file="figures/1d-cov-structures.pdf", width=4, height=2.5)
par(mar=c(3,3,1,1))
plot_samples(x_seq, se_samples, ylim=c(-2,2), lwd=0.5)
lines(x_seq, as.vector(se_samples[2,]), col="red3", lwd=2)

plot_samples(x_seq, matern_samples, ylim=c(-2,2), lwd=0.5)
lines(x_seq, as.vector(matern_samples[5,]), col="red3", lwd=2)

plot_samples(x_seq, exp_samples, ylim=c(-2,2), lwd=0.5)
lines(x_seq, as.vector(exp_samples[3,]), col="red3", lwd=1.2)
dev.off()


# =============================================
# === 2D PLOTS
# =============================================

n <- 81
x_seq <- seq(0, 1, l=n)
x_lat <- expand.grid(x1=x_seq, x2=x_seq)

dmat <- as.matrix(dist(x_lat/theta))

se_cov_mat <- sigma^2 * se_kernel(dmat)
matern_cov_mat <- sigma^2 * matern_kernel(dmat)
exp_cov_mat <- sigma^2 * exp_kernel(dmat)

m <- rep(0, n^2)

set.seed(42)
# Caution: these may take a while to run, optionally can reduce `n` to speed this up
se_samples <- rmvnorm(n=1, mean=m, sigma=se_cov_mat)
matern_samples <- rmvnorm(n=1, mean=m, sigma=matern_cov_mat)
exp_samples <- rmvnorm(n=1, mean=m, sigma=exp_cov_mat)

z_se <- matrix(se_samples[1,], nrow=n, ncol=n)
z_matern <- matrix(matern_samples[1,], nrow=n, ncol=n)
z_exp <- matrix(exp_samples[1,], nrow=n, ncol=n)

pdf(file="figures/2d-cov-structures.pdf", width=4, height=4)
par(mar=c(3,3,1,1))
#filled.contour(x_seq, x_seq, z_se, color.palette=function(n){hcl.colors(n, "YlOrRd")}, nlevels=20)
cpal <- function(n) hcl.colors(n, "YlOrRd")
image(x_seq, x_seq, z_se, col=cpal(20), axes=T, ylab="", xlab="", xlim=c(0,1), ylim=c(0,1)); box()
image(x_seq, x_seq, z_matern, col=cpal(20), axes=T, ylab="", xlab="", xlim=c(0,1), ylim=c(0,1)); box()
image(x_seq, x_seq, z_exp, col=cpal(100), axes=T, ylab="", xlab="", xlim=c(0,1), ylim=c(0,1)); box()
dev.off()
