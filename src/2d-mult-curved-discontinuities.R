library(viridisLite)
library(plot3D)

### Plotting setup
axis_labels <- c("", "", "f(x,y)")  #c(expression(x), expression(y), expression(f(x,y)))
points_plot <- function(xD, col="black", cex=1.5, pch=20, ...) points(xD, col=col, cex=cex, pch=pch, ...)

var_cols <- function(n) hcl.colors(n, "Temps")  #function(n) hcl.colors(n, "YlOrRd", rev=TRUE)
surf_cols <- plasma  #magma
diag_cols <- turbo

discontinuity_plot <- function(col="black", lwd=4) {
  xseq <- seq(-1, 1, length=500)
  yseq1 <- -(0.2*exp((-xseq-0.8)^2) - 0.85)
  ind1 <- which(yseq1 > -0.25)
  x1 <- xseq[ind1]; y1 <- yseq1[ind1]
  lines(x1, y1, col=col, lwd=lwd)
  
  yseq2 <- (0.2*exp((xseq-0.8)^2) - 0.85)
  ind2 <- which(yseq2 < 0.25)
  x2 <- xseq[ind2]; y2 <- yseq2[ind2]
  lines(x2, y2, col=col, lwd=lwd)
}


### ===================================
### Embedding Surface
### ===================================


# Embedding surface
v <- function(X){
  if (is.matrix(X)) { x <- X[,1]; y <- X[,2] }
  else { x <- X[1]; y <- X[2] }
  
  A1 <- 0.2*exp((x-0.8)^2) - 0.85
  A2 <- 0.2*exp((-x-0.8)^2) - 0.85
  x0 <- -0.50566
  
  v1 <- -(y < 0.25)*(y < A1)*( 0.5*(x > x0)*(x - x0)^2 + (y - 0.25)^2 )
  v2 <- (y > -0.25)*(y > -A2)*( 0.5*(x < -x0)*(x + x0)^2 + (y + 0.25)^2 )
  
  return(v1 + v2)
}


V <- function(X) cbind(X, v(X))


# Plot embedding surface
n_seq <- 501; x_seq <- seq(-1,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))
x_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=FALSE)
y_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=TRUE)

vP <- v(xP)
vP_mat <- matrix(vP, nrow=n_seq, ncol=n_seq)
vP_mat[which((y_mat > -0.25) & (abs(y_mat + (0.2*exp((-x_mat-0.8)^2) - 0.85)) < 0.005))] <- NA
vP_mat[which((y_mat < 0.25) & (abs(y_mat - (0.2*exp((x_mat-0.8)^2) - 0.85)) < 0.005))] <- NA


#filled.contour(x=x_seq, y=x_seq, z=vP_mat, color.palette=magma, levels=seq(-3, 3, 0.1),
#               xlab=axis_labels[1], ylab=axis_labels[2],
#               plot.axes={axis(1);axis(2)
#                 discontinuity_plot()
#               })


png(filename="figures/disc2-embedding-surf-3d.png", width=6, height=6, units='in', res=600, bg="transparent")
par(mar=c(2,2,2,2))
persp3D(x_seq, x_seq, vP_mat,
        theta=20, phi=20, expand=0.8, colkey=FALSE, col=magma(200)[10:190], shade=2, facets=F,
        ltheta=20, lphi=50, zlab="v(x,y)",
        lighting=list(ambient=0.4, diffuse=0.6, specular=0.3, exponent=20, sr=0, alpha=1), border=NA, lwd=2)
dev.off()




### True Function
f <- function(X){
  if (is.matrix(X)) {
    x <- X[,1]; y <- X[,2]
  } else {
    x <- X[1]; y <- X[2]
  }
  
  V <- v(X)
  A <- V*sin(x*y - 4)^2 *y + 5*(sqrt(cos(x^2 - y^2)) - 1.1)
  E <- exp(-5*((x+0.5)^2 + (y+0.3)^2))*(2*(y < (0.2*exp((x-0.8)^2) - 0.85)) - 1)*(y-0.25)*(y<0.25)
  B <- tanh(10*exp(-3*(x^4 + y^4)))
  
  return((2*tanh(abs(A + 0.8)) - 0.8 + E)*B)
}



n_seq <- 501
x_seq <- seq(-1,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))

# Evaluate true function
fP <- f(xP)
fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq) # resize to matrix for plotting


### True function plot
# 2D contour plot
pdf(file="figures/disc2-true-2d.pdf", width=7, height=6)
par(mar=c(4,4,2,4))
filled.contour(x=x_seq, y=x_seq, z=fP_mat, color.palette=surf_cols, levels=seq(-1.5, 1.5, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, fP_mat, add=TRUE, level=seq(-1.5, 1.5, 0.1), lwd=0.4, drawlabels=FALSE)
                 discontinuity_plot(lwd=5)
               })
dev.off()


# 3D plot
n_seq <- 501; x_seq <- seq(-1,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))
fP <- f(xP); fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)
x_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=FALSE)
y_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=TRUE)

# Mask Discontinuity
fP_mat[which((y_mat > -0.25) & (abs(y_mat + (0.2*exp((-x_mat-0.8)^2) - 0.85)) < 0.005))] <- NA
fP_mat[which((y_mat < 0.25) & (abs(y_mat - (0.2*exp((x_mat-0.8)^2) - 0.85)) < 0.005))] <- NA

png(filename="figures/disc2-true-3d.png", width=6, height=6, units='in', res=600, bg="transparent")
par(mar=c(2,2,2,2))
persp3D(x_seq, x_seq, fP_mat,
        theta=15, phi=35, expand=0.6, colkey=FALSE, col=plasma(100), shade=1, facets=FALSE,
        ltheta=20, lphi=50, zlab=axis_labels[3],
        lighting=list(ambient=0.5, diffuse=0.6, specular=0.3, exponent=20, sr=0, alpha=1), border=NA, lwd=2)
dev.off()





# Design
nD_seq <- 6
xD_seq <- seq(-1+1/nD_seq, 1-1/nD_seq, len=nD_seq)

xD <- as.matrix(expand.grid(xD_seq, xD_seq))
fD <- f(xD)

n_seq <- 201; x_seq <- seq(-1,1, len=n_seq) # may be slow to run emulator and can decrease this to 101
xP <- as.matrix(expand.grid(x_seq, x_seq))
x_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=FALSE)
y_mat <- matrix(x_seq, nrow=n_seq, ncol=n_seq, byrow=TRUE)

fP <- f(xP)
fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)


# Embed in Higher Dimension
xP_h <- V(xP)
xD_h <- V(xD)



### ===================================
### TENSE
### ===================================


v_deriv <- function(X){
  if (is.matrix(X)) { x <- X[,1]; y <- X[,2] }
  else { x <- X[1]; y <- X[2] }
  
  A1 <- 0.2*exp((x-0.8)^2) - 0.85
  A2 <- 0.2*exp((-x-0.8)^2) - 0.85
  x0 <- 0.50566
  
  v1_x_deriv <- -(y < 0.25)*(y < A1)*(x > -x0)*(x + x0)
  v1_y_deriv <- -(y < 0.25)*(y < A1)*2*(y - 0.25)
  
  v2_x_deriv <- (y > -0.25)*(y > -A2)*(x < x0)*(x - x0)
  v2_y_deriv <- (y > -0.25)*(y > -A2)*2*(y + 0.25)
  
  x_deriv <- v1_x_deriv + v2_x_deriv
  y_deriv <- v1_y_deriv + v2_y_deriv
  
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

### Matern 3/2 covariance function
k_S   <- function(d, sig=1, nu=5/2, theta=1) {
  A = 2*sqrt(nu) * d/theta^2
  return(sig^2 * 2^(1-nu) / gamma(nu) * (A)^nu * besselK(A, nu))
}
k_NS  <- function(x, y, sig=1, Sigma) {  # non-stationary covariance function
  p <- length(x)
  2^(p/2) * det(Sigma(x))^(1/4) * det(Sigma(y))^(1/4) / det(Sigma(x)+Sigma(y))^(1/2) * k_S(d=sqrt(Q(x,y)), sig=sig)
}


source("src/lib/simple_NS_emulator.R")



#####========================================
### TENSE Emulation
#####========================================

NS_em <- simple_NS_emulator(xD_h, fD, xP_h, sig=1, nugget=1e-3, mu=0, Sigma=Sigma_real, just_var=TRUE)
exp <- NS_em[["ExpD_f(x)"]]
var <- NS_em[["VarD_f(x)"]]

exp_mat <- matrix(exp, nrow=n_seq, ncol=n_seq)
var_mat <- matrix(var, nrow=n_seq, ncol=n_seq)
sd_mat <- sqrt(abs(var_mat))
diag_mat <- (exp_mat - fP_mat)/sd_mat

# Emulator expectation plot
pdf(file="figures/disc2-TENSE-2d.pdf", width=7, height=6)
par(mar=c(4,4,2,4))
filled.contour(x=x_seq, y=x_seq, z=exp_mat, color.palette=surf_cols, level=seq(-1.5, 1.5, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, exp_mat, add=TRUE, level=seq(-1.5, 1.5, 0.1), lwd=0.4, drawlabels=FALSE)
                 discontinuity_plot()
                 points_plot(xD)
               })

# Emulator sd plot
filled.contour(x=x_seq, y=x_seq, z=sd_mat, color.palette=var_cols, levels=seq(0,1,0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 discontinuity_plot()
                 points_plot(xD)
               })

# Emulator Diagnostics
filled.contour(x=x_seq, y=x_seq, z=diag_mat, color.palette=diag_cols, levels=seq(-3,3,0.25),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, diag_mat, add=TRUE, level=seq(-3,3,0.25), lwd=0.4, drawlabels=FALSE)
                 discontinuity_plot()
                 points_plot(xD)
               })
dev.off()

