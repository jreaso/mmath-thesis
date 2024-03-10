library('SimDesign')

### GP PLOTTING FUNCTION
# Reuses code from `1d-gp.R`, but modified

plot_gp <- function(xP, f, samples=NULL, ylim=NULL){
  int_col <- rgb(.9,.9,.9)
  sample_col <- rgb(.7,.7,.7)
  
  # Initialise the plot window
  plot(1, 1, xlim=c(0,1), ylim=ylim, type="n", xlab="", ylab="")
  
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

# Simple Non-Stationary Squared Exponential Kernel
theta <- function(x) {
  0.2 - exp(-200*(x-0.5)^2) * 3/16
}

cf <- function(x, y) {
  sqrt((2*theta(x)*theta(y))/(theta(x)^2 + theta(y)^2)) * exp(-(x-y)^2/(theta(x)^2 + theta(y)^2))
}

xP <- seq(-.05, 1.05, len=1001)

gp_exp_prior <- rep(0, length(xP))
gp_var_prior_mat <- outer(xP, xP, cf)

# Sample from posterior
set.seed(4)
n_samples <- 3
prior_samples <- rmvnorm(n=n_samples, mean=gp_exp_prior, sigma=gp_var_prior_mat)

pdf(file="figures/1d-simple-nscs.pdf", width=6, height=5)
  par(mar = c(4, 4, 2, 2))
  plot_gp(f=theta(xP), xP=xP, ylim=c(0, 0.23))
  title(ylab=expression(theta(x)), xlab="x", line=2.5, cex.lab=1)
  
  
  plot_gp(f=prior_samples[1,], xP=xP, 
          samples=prior_samples[2:3,], ylim=c(-2.5, 2.5))
  title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)
dev.off()

