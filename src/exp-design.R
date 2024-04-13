library(LHD)
library(lhs)

# 2D maximin LHD - we mostly use library for efficiency
source("src/lib/maximin_lhd_2d.R")


### Function to plot 2D design
plot_design <- function(X, grid=T, grid_col="grey60", xlim=c(0,1), ylim=c(0,1), pch=16,
                        xaxs="i", yaxs="i", col="dodgerblue3", xlab="", ylab="", cex=1.2){
  plot(X, xlim=xlim, ylim=ylim, pch=pch, xaxs=xaxs, yaxs=yaxs, col=col, xlab=xlab, ylab=ylab, cex=cex)
  if (grid){
    nl <- dim(X)[1]
    abline(h=(0:nl)/nl, col=grid_col)
    abline(v=(0:nl)/nl, col=grid_col)
  }
}


pdf(file="figures/2d-random-lhd.pdf", width=4, height=4)
par(mar=c(4,4,2,2))
# Random LHD
set.seed(5)
X_r <- rLHD(n=5, k=2)/5 - 1/10
plot_design(X_r, cex=2)

# Random Perturbed LHD
X_r_p <- X_r + matrix(rnorm(10, mean=0, sd=1/20), nrow=5, ncol=2)
plot_design(X_r_p, cex=2)
dev.off()



pdf(file="figures/2d-Mm-AE-lhd.pdf", width=5, height=5)
par(mar=c(4,4,2,2))
# Maximin LHD
X_Mm1 <- maximin_lhd_2d(16, perturb=F)
plot_design(X_Mm1)

#X_Mm2 <- FastMmLHD(n=16, k=2)/16 - 1/32
#plot_design(X_Mm2)

#X_Mm3 <- maximinLHS(n=16, k=2, method="iterative", dup=2, maxIter=1000, eps=1e-12, optimize.on="result", debug=T)
#plot_design(X_Mm3)
  

# Genetic Algorithm LHD
#X_g <- GA(n=16, k=2)/16 - 1/32
#plot_design(X_g)



# Audze-Eglaise LHD - from space filling designs
# d = 2.6659008
X_ae <- matrix(c(0, 10, 1, 3, 2, 14, 3, 7, 4, 0, 5, 11, 6, 4, 7, 15, 8, 8, 9, 1, 10, 12, 11, 5, 12, 9, 13, 2, 14, 13, 15, 6), nrow=16, ncol=2, byrow=T)/16 - 1/32
plot_design(X_ae)
#X_ae
dev.off()

