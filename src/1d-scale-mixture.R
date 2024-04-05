library('SimDesign')

### GP PLOTTING FUNCTION
# Reuses code from `1d-gp.R`, but modified

plot_gp <- function(xP, f, samples=NULL, ylim=NULL, xlim=c(0,1)){
  int_col <- rgb(.9,.9,.9)
  sample_col <- rgb(.7,.7,.7)
  
  # Initialise the plot window
  plot(1, 1, xlim=xlim, ylim=ylim, type="n", xlab="", ylab="")
  
  # Plot samples
  if (!is.null(samples)) {
    for (i in 1:dim(samples)[1]){
      lines(xP, as.vector(samples[i,]), col=sample_col, lwd=1)
    }
  }
  
  # Plot the highltighted function
  lines(xP, f, lwd=2.5)
  
  box()
}


### PRIOR PLOT

# Scale mixture of gaussians
k <- function(d) {
  d2 <- d^2
  return(0.25*exp(-d2/0.08^2) + exp(-d2/0.8^2))
}

cf <- function(x, y) k((x-y))



xP <- seq(-.05, 1.05, len=1001)
x_seq <- seq(-1.05, 1.05, len=1001)

gp_exp_prior <- rep(0, length(xP))
gp_var_prior_mat <- outer(xP, xP, cf)

# Sample from posterior
set.seed(1)
n_samples <- 3
prior_samples <- rmvnorm(n=n_samples, mean=gp_exp_prior, sigma=gp_var_prior_mat)

pdf(file="figures/1d-scale-mixture.pdf", width=6, height=4)
par(mar = c(4, 4, 2, 2))
plot_gp(f=k(x_seq), xP=x_seq, ylim=c(0, 1.5), xlim=c(-1,1))
title(ylab="k(x, x')", xlab="x-x'", line=2.5, cex.lab=1)


plot_gp(f=prior_samples[1,], xP=xP, 
        samples=prior_samples[2:3,], ylim=c(-2.5, 2.5))
title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)
dev.off()

