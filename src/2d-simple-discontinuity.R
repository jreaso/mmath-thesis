library(viridisLite)
library(plot3D)

# Simple BL Emulator
source("src/lib/simple_BL_emulator.R")


### Plotting setup
axis_labels <- c("", "", "f(x,y)")  #c(expression(x), expression(y), expression(f(x,y)))
points_plot <- function(xD, col="black", cex=1.5, pch=20, ...) points(xD, col=col, cex=cex, pch=pch, ...)

var_cols <- function(n) hcl.colors(n, "Temps")  #function(n) hcl.colors(n, "YlOrRd", rev=TRUE)
surf_cols <- plasma  #magma
diag_cols <- turbo

discontinuity_plot <- function(col="black", lwd=4) {
  lines(c(0.5, 0.5), c(0, 0.5), col=col, lwd=lwd)
}


### True Function
f <- function(X){
  if (is.matrix(X)) {
    x <- X[,1]; y <- X[,2]
  } else {
    x <- X[1]; y <- X[2]
  }
  
  t <- (y < 0.5)*(y - 0.5)^2 * (1 - 2*(x < 0.5))
  A1 <- exp(-10*((x-0.4)^2 + (y-0.7)^2))
  A2 <- exp(-20*((x-0.9)^2 + (y-0.5)^2))
  B <- cos(8*((x-0.5)^2 - y^2))
  
  Z <- t*5*(sin(3*x)^2 + x) + B/2 + 2*A1 + 2*A2
  
  return(-Z/3 + 0.35)
}



n_seq <- 101
x_seq <- seq(0,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))

# Evaluate true function
fP <- f(xP)
fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq) # resize to matrix for plotting


### True function plot
# 2D contour plot
pdf(file="figures/disc1-true-2d.pdf", width=7, height=6)
par(mar-c(4,4,2,4))
filled.contour(x=x_seq, y=x_seq, z=fP_mat, color.palette=surf_cols, levels=seq(-1, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, fP_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                 discontinuity_plot()
               })
dev.off()

# 2D contour plot with geodesic lines
pdf(file="figures/disc1-true-2d-geodesic.pdf", width=7, height=6)
par(mar-c(4,4,2,4))
filled.contour(x=x_seq, y=x_seq, z=fP_mat, color.palette=surf_cols, levels=seq(-1, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, fP_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                 discontinuity_plot()
                 
                 # Geodesic
                 lines(x=c(0.4, 0.5, 0.65), y=c(0.05, 0.5, 0.2), col='white', lwd=3)
                 lines(x=c(0.4, 0.65), y=c(0.05, 0.2), col='white', lwd=1.5, lty=2)
                 points(x=c(0.4, 0.5, 0.65), y=c(0.05, 0.5, 0.2),
                        col="black", pch=21, bg="black", cex=1.4)
               })
dev.off()

# 3D plot
n_seq <- 501; x_seq <- seq(0,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))
fP <- f(xP); fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)
x_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=FALSE)
y_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=TRUE)

# Mask Discontinuity
fP_mat[which(x_mat > 0.498 & x_mat < 0.502 & y_mat < 0.5)] <- NA

png(filename="figures/disc1-true-3d.png", width=6, height=6, units='in', res=600, bg="transparent")
  par(mar=c(2,2,2,2))
  persp3D(x_seq, x_seq, fP_mat,
          theta=25, phi=25, expand=0.6, colkey=FALSE, col=plasma(100), shade=1, facets=FALSE,
          ltheta=20, lphi=50, zlab=axis_labels[3],
          lighting=list(ambient=0.5, diffuse=0.6, specular=0.3, exponent=20, sr=0, alpha=1), border=NA, lwd=2)
dev.off()



### Standard Emulator (smooths over discontinuity)

n_seq <- 501; x_seq <- seq(0,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))

fP <- f(xP); fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)

n1_seq <- 32
x1_seq <- seq(0+1/(2*n1_seq), 1-1/(2*n1_seq), len=n1_seq)

xD1 <- as.matrix(expand.grid(x1_seq, x1_seq))
fD1 <- f(xD1)

# Warning: this may take a while to run
BL_em1 <- simple_BL_emulator(xD=xD1, D=fD1, xP=xP, theta=0.6, sig=5, mu=0, nugget=1e-2)  
exp1 <- BL_em1[["ExpD_f(x)"]]
var1 <- BL_em1[["VarD_f(x)"]]
exp1_mat <- matrix(exp1, nrow=n_seq, ncol=n_seq)
var1_mat <- matrix(var1, nrow=n_seq, ncol=n_seq)

png(filename="figures/disc1-standard-smooth-3d.png", width=6, height=6, units='in', res=600, bg="transparent")
  par(mar=c(2,2,2,2))
  persp3D(x_seq, x_seq, exp1_mat,
          theta=25, phi=25, expand=0.6, colkey=FALSE, col=plasma(100), shade=1, facets=F,
          ltheta=20, lphi=50, zlab=axis_labels[3],
          lighting=list(ambient=0.5, diffuse=0.6, specular=0.3, exponent=20, sr=0, alpha=1), border=NA, lwd=2)
dev.off()


### INCOMPLETE - need to perform torn embedding emulation and TENSE emulation
