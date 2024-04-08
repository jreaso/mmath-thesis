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


# Random LHD
X_r <- rLHD(n=16, k=2)/16 - 1/32
plot_design(X_r)


# Maximin LHD
X_Mm1 <- maximin_lhd_2d(16, perturb=F)
plot_design(X_Mm1)

#X_Mm2 <- FastMmLHD(n=16, k=2)/16 - 1/32
#plot_design(X_Mm2)

X_Mm3 <- maximinLHS(n=16, k=2, method="iterative", maxIter=1000, eps=1e-10, optimize.on="result", debug=T)
plot_design(X_Mm3)
  

# Genetic Algorithm LHD
X_g <- GA(n=16, k=2)/16 - 1/32
plot_design(X_g)

