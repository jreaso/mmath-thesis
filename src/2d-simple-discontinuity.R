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
#pdf(file="figures/disc1-true-2d-geodesic.pdf", width=7, height=6)
#par(mar-c(4,4,2,4))
#filled.contour(x=x_seq, y=x_seq, z=fP_mat, color.palette=surf_cols, levels=seq(-1, 1, 0.1),
#               xlab=axis_labels[1], ylab=axis_labels[2],
#               plot.axes={axis(1);axis(2)
#                 contour(x_seq, x_seq, fP_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
#                 discontinuity_plot()
#                 
#                 # Geodesic
#                 lines(x=c(0.4, 0.5, 0.65), y=c(0.05, 0.5, 0.2), col='white', lwd=5)
#                 lines(x=c(0.4, 0.65), y=c(0.05, 0.2), col='white', lwd=4, lty=2)
#                 points(x=c(0.4, 0.5, 0.65), y=c(0.05, 0.5, 0.2),
#                        col="black", pch=21, bg="black", cex=2)
#               })
#dev.off()

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


### Pseudo NS CS
# Instead of building a NS CS, we build a function slightly which is f with a smoothed over discontinuity by using a sigmoid

f_smooth <- function(X){
  if (is.matrix(X)) {
    x <- X[,1]; y <- X[,2]
  } else {
    x <- X[1]; y <- X[2]
  }
  
  t <- (y < 0.5)*(y - 0.5)^2 * (2/(1 + exp(-100*(x-0.5))) - 1)
  A1 <- exp(-10*((x-0.4)^2 + (y-0.7)^2))
  A2 <- exp(-20*((x-0.9)^2 + (y-0.5)^2))
  B <- cos(8*((x-0.5)^2 - y^2))
  
  Z <- t*5*(sin(3*x)^2 + x) + B/2 + 2*A1 + 2*A2
  
  return(-Z/3 + 0.35)
}

# 3D plot
n_seq <- 1001; x_seq <- seq(0,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))
fsP <- f_smooth(xP); fsP_mat <- matrix(fsP, nrow=n_seq, ncol=n_seq)
x_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=FALSE)
y_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=TRUE)

png(filename="figures/disc1-pseudo-nscs-3d.png", width=6, height=6, units='in', res=600, bg="transparent")
par(mar=c(2,2,2,2))
persp3D(x_seq, x_seq, fsP_mat,
        theta=25, phi=25, expand=0.6, colkey=FALSE, col=plasma(100), shade=1, facets=F,
        ltheta=20, lphi=50, zlab=axis_labels[3],
        lighting=list(ambient=0.5, diffuse=0.6, specular=0.3, exponent=20, sr=0, alpha=1), border=NA, lwd=2)
dev.off()



### ===================================
### Stationary Torn Embedding Emulation
### ===================================


# Embedding surface
v <- function(X){
  if (is.matrix(X)) {
    x <- X[,1]; y <- X[,2]
  } else {
    x <- X[1]; y <- X[2]
  }
  (y < 0.5)*(y - 0.5)^2 * (2*(x < 0.5) - 1)
}

V <- function(X) cbind(X, v(X))


# Plot embedding surface
n_seq <- 501; x_seq <- seq(0,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))
x_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=FALSE)
y_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=TRUE)

vP <- v(xP)
vP_mat <- matrix(vP, nrow=n_seq, ncol=n_seq)
vP_mat[which(x_mat > 0.498 & x_mat < 0.502 & y_mat < 0.5)] <- NA


png(filename="figures/disc1-embedding-surf-3d.png", width=6, height=6, units='in', res=600, bg="transparent")
par(mar=c(2,2,2,2))
persp3D(x_seq, x_seq, vP_mat,
        theta=40, phi=20, expand=0.6, colkey=FALSE, col=magma(100), shade=1, facets=F,
        ltheta=20, lphi=50, zlab="v(x,y)",
        lighting=list(ambient=0.5, diffuse=0.6, specular=0.3, exponent=20, sr=0, alpha=1), border=NA, lwd=2)
dev.off()


# Design
n2_seq <- 4
x2_seq <- seq(0+1/(2*n2_seq), 1-1/(2*n2_seq), len=n2_seq)

xD2 <- as.matrix(expand.grid(x2_seq, x2_seq))
fD2 <- f(xD2)

n_seq <- 101; x_seq <- seq(0,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))
x_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=FALSE)
y_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=TRUE)


# Embed in Higher Dimension
xP_h <- V(xP)
xD2_h <- V(xD2)


# Torn Embedding Stationary Emulation
BL_em2 <- simple_BL_emulator(xD2_h, fD2, xP_h, theta=0.3, sig=1, nugget=1e-3, mu=0)
exp2 <- BL_em2[["ExpD_f(x)"]]
var2 <- BL_em2[["VarD_f(x)"]]

exp2_mat <- matrix(exp2, nrow=n_seq, ncol=n_seq)
var2_mat <- matrix(var2, nrow=n_seq, ncol=n_seq)
sd2_mat <- sqrt(abs(var2_mat))  # absolute taken for stability

# Emulator expectation plot
pdf(file="figures/disc1-stat-torn-embed.pdf", width=7, height=6)
par(mar=c(4,4,2,4))
  filled.contour(x=x_seq, y=x_seq, z=exp2_mat, color.palette=surf_cols, level=seq(-1, 1, 0.1),
                 xlab=axis_labels[1], ylab=axis_labels[2],
                 plot.axes={axis(1);axis(2)
                   contour(x_seq, x_seq, exp2_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                   discontinuity_plot()
                   points_plot(xD2)
                 })

# Emulator sd plot
  filled.contour(x=x_seq, y=x_seq, z=sd2_mat, color.palette=var_cols, levels=seq(0,1,0.1),
                 xlab=axis_labels[1], ylab=axis_labels[2],
                 plot.axes={axis(1);axis(2)
                   discontinuity_plot()
                   points_plot(xD2)
                 })


# Use a more extreme embedding surface to highlight warped surfaces
v_extreme <- function(X){
  if (is.matrix(X)) {
    x <- X[,1]; y <- X[,2]
  } else {
    x <- X[1]; y <- X[2]
  }
  3*(y < 0.5)*(y - 0.5)^2 * (2*(x < 0.5) - 1)
}

V_extreme <- function(X) cbind(X, v_extreme(X))

xP_h <- V_extreme(xP)
xD2_h <- V_extreme(xD2)


# Torn Embedding Stationary Emulation
BL_em3 <- simple_BL_emulator(xD2_h, fD2, xP_h, theta=0.3, sig=1, nugget=1e-3, mu=0)
exp3 <- BL_em3[["ExpD_f(x)"]]
var3 <- BL_em3[["VarD_f(x)"]]

exp3_mat <- matrix(exp3, nrow=n_seq, ncol=n_seq)
var3_mat <- matrix(var3, nrow=n_seq, ncol=n_seq)
sd3_mat <- sqrt(abs(var3_mat))  # absolute taken for stability

# Emulator expectation plot
  filled.contour(x=x_seq, y=x_seq, z=exp3_mat, color.palette=surf_cols, level=seq(-1, 1, 0.1),
                 xlab=axis_labels[1], ylab=axis_labels[2],
                 plot.axes={axis(1);axis(2)
                   contour(x_seq, x_seq, exp3_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                   discontinuity_plot()
                   points_plot(xD2)
                 })

# Emulator sd plot
  filled.contour(x=x_seq, y=x_seq, z=sd3_mat, color.palette=var_cols, levels=seq(0,1,0.1),
                 xlab=axis_labels[1], ylab=axis_labels[2],
                 plot.axes={axis(1);axis(2)
                   discontinuity_plot()
                   points_plot(xD2)
                 })
dev.off()



### ===================================
### TENSE
### ===================================

# Work with v_extreme embedding surface

v_deriv <- function(X){
  if (is.matrix(X)) {
    x <- X[,1]; y <- X[,2]
  } else {
    x <- X[1]; y <- X[2]
  }
  x_deriv <- 0
  y_deriv <- 6*(y < 0.5)*(y - 0.5) * (2*(x < 0.5) - 1)
  
  return(c(x_deriv, y_deriv))
}


Sigma_real <- function(X, theta=0.3, alpha_3=0.25) {
  v_d <- v_deriv(X)
  v_x <- v_d[1]; v_y <- v_d[2]
  w_3 <- c(-v_x, -v_y, 1)
  r2 <- v_x^2 + v_y^2
  return(										
    theta^2 * matrix(c(  1,   0,   v_x,
                         0,   1,   v_y,
                         v_x, v_y,    r2), nrow=3, ncol=3, byrow=TRUE) + alpha_3^2 / (r2+1) * w_3 %*% t(w_3)
  )
}

### TENSE NS CS
Q 	  <- function(x, y, Sigma=Sigma_real) (x-y) %*% solve( (Sigma(x) + Sigma(y))/2 ) %*% (x-y)  # Quadratic Form
k_S   <- function(d, sig=1) sig^2*exp(-d^2)  # squared exponential covariance function
k_NS  <- function(x, y, sig=1, Sigma) {  # non-stationary covariance function
  p <- length(x)
  2^(p/2) * det(Sigma(x))^(1/4) * det(Sigma(y))^(1/4) / det(Sigma(x)+Sigma(y))^(1/2) * k_S(d=sqrt(Q(x,y)), sig=sig)
  
}


source("src/lib/simple_NS_emulator.R")



#####========================================
### TENSE NS Emulation
#####========================================

NS_em4 <- simple_NS_emulator(xD2_h, fD2, xP_h, sig=1, nugget=1e-3, mu=0, Sigma=Sigma_real, just_var=TRUE)
exp4 <- NS_em4[["ExpD_f(x)"]]
var4 <- NS_em4[["VarD_f(x)"]]

exp4_mat <- matrix(exp4, nrow=n_seq, ncol=n_seq)
var4_mat <- matrix(var4, nrow=n_seq, ncol=n_seq)
sd4_mat <- sqrt(abs(var4_mat))

# Emulator expectation plot
pdf(file="figures/disc1-TENSE-2d.pdf", width=7, height=6)
par(mar=c(4,4,2,4))
  filled.contour(x=x_seq, y=x_seq, z=exp4_mat, color.palette=surf_cols, level=seq(-1, 1, 0.1),
                 xlab=axis_labels[1], ylab=axis_labels[2],
                 plot.axes={axis(1);axis(2)
                   contour(x_seq, x_seq, exp4_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                   discontinuity_plot()
                   points_plot(xD2)
                 })

# Emulator sd plot
filled.contour(x=x_seq, y=x_seq, z=sd4_mat, color.palette=var_cols, levels=seq(0,1,0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 discontinuity_plot()
                 points_plot(xD2)
               })
dev.off()

