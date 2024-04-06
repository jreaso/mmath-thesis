library(SimDesign)


# =============================================
# === COVARIANCE STRUCTURES
# =============================================

sigma <- 1; theta <- 0.35; p <- 0.45

# Periodic kernel
periodic_kernel <- function(dist_mat) exp(-(sin(dist_mat * pi/p))^2 / theta^2)

# Locally periodic kernel
locally_periodic_kernel <- function(dist_mat) exp(-(sin(dist_mat * pi/p))^2 / theta^2) * exp(-(dist_mat)^2 / theta^2)


# =============================================
# === 1D PLOTS
# =============================================

n <- 1001
x_seq <- seq(-0.1, 1.1, l=n)

dmat <- as.matrix(dist(x_seq))

p_cov_mat <- sigma^2 * periodic_kernel(dmat)
lp_cov_mat <- sigma^2 * locally_periodic_kernel(dmat)

set.seed(1)
p_samples <- rmvnorm(n=4, mean=rep(0, n), sigma=p_cov_mat)
lp_samples <- rmvnorm(n=4, mean=rep(0, n), sigma=lp_cov_mat)


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


pdf(file="figures/1d-periodic-cf.pdf", width=6, height=4)
par(mar=c(3,3,1,1))
plot_samples(x_seq, p_samples, ylim=c(-2,2), lwd=0.5)
lines(x_seq, as.vector(p_samples[1,]), col="red3", lwd=2)

plot_samples(x_seq, lp_samples, ylim=c(-2,2), lwd=0.5)
lines(x_seq, as.vector(lp_samples[2,]), col="red3", lwd=2)
dev.off()
