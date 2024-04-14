library(viridisLite)
library(plotly)

source("src/lib/filled.contour3.R")
source("src/lib/TENSE.R")

###====================================
### Embedding Surface
###====================================

# Compact domain sigmoid function on [0, 1]
s <- function(d) {
  (2*d^2 * (0 <= d & d < 0.5)) + ((1 - 2*(d-1)^2) * (0.5 <= d & d < 1)) + (1*(1 <= d))
}

# Embedding surface
v <- function(X){
  if (is.matrix(X)) { x <- X[,1]; y <- X[,2]; z <- X[,3] }
  else { x <- X[1]; y <- X[2]; z <- X[3] }
  
  A <- 1 - s(3*sqrt((x-0.5)^2 + (0.5*z-0)^2 + 1e-6))
  B <- 1 - 2*(y > 0.5)
  return(A*B)
}

# Embedding function
V <- function(X) cbind(X, v(X))



###====================================
### Plotting Embedding Surface
###====================================
n_seq <- 201
x_seq <- seq(0,1, l=n_seq);
x2_grid <- as.matrix(expand.grid(x_seq, x_seq))

#pdf(file="figures/T3D-....pdf", width=7, height=6)
xP1 <- cbind(x2_grid[,1], 0.25, x2_grid[,2]); vP1 <- v(xP1)
vP1_mat <- matrix(vP1, nrow=n_seq, ncol=n_seq)
par(mar=c(4,4,2,2))
filled.contour3(x=x_seq, y=x_seq, z=vP1_mat, color.palette=magma, levels=seq(-1, 1, 0.05), xlab="x", ylab="z")

xP2 <- cbind(x2_grid[,1], 0.75, x2_grid[,2]); vP2 <- v(xP2)
vP2_mat <- matrix(vP2, nrow=n_seq, ncol=n_seq)
filled.contour3(x=x_seq, y=x_seq, z=vP2_mat, color.palette=magma, levels=seq(-1, 1, 0.05), xlab="x", ylab="z")

xP3 <- cbind(0.5, x2_grid[,1], x2_grid[,2]); vP3 <- v(xP3)
vP3_mat <- matrix(vP3, nrow=n_seq, ncol=n_seq)
filled.contour3(x=x_seq, y=x_seq, z=vP3_mat, color.palette=magma, levels=seq(-1, 1, 0.05),
                xlab="y", ylab="z",
                plot.axes={axis(1);axis(2)
                  lines(c(0.5, 0.5), c(0, 2/3), type="l", lwd=3)
                })

xP4 <- cbind(0.35, x2_grid[,1], x2_grid[,2]); vP4 <- v(xP4)
vP4_mat <- matrix(vP4, nrow=n_seq, ncol=n_seq)
filled.contour3(x=x_seq, y=x_seq, z=vP4_mat, color.palette=magma, levels=seq(-1, 1, 0.05),
                xlab="y", ylab="z",
                plot.axes={axis(1);axis(2)
                  lines(c(0.5, 0.5), c(0, 0.6), type="l", lwd=3)
                })

xP5 <- cbind(0.3, x2_grid[,1], x2_grid[,2]); vP5 <- v(xP5)
vP5_mat <- matrix(vP5, nrow=n_seq, ncol=n_seq)
filled.contour3(x=x_seq, y=x_seq, z=vP5_mat, color.palette=magma, levels=seq(-1, 1, 0.05),
                xlab="y", ylab="z",
                plot.axes={axis(1);axis(2)
                  lines(c(0.5, 0.5), c(0, 0.533), type="l", lwd=3)
                })
#dev.off()



###====================================
### True Function
###====================================

f <- function(X){
  if (is.matrix(X)) { x <- X[,1]; y <- X[,2]; z <- X[,3] }
  else { x <- X[1]; y <- X[2]; z <- X[3] }
  
  vX <- v(X)
  
  # basis: 1, x, y, z, x^2, y^2, z^2, xz
  R <- 3*(x-0.4)^2 + 2*y^2 + 5*(x-0.3)*(z-0.3) - 3.5
  
  G1 <- exp(-50*((x-0.5)^2 + (z-0.6)^2))
  G2 <- exp(-30*((x-0.8)^2 + (z-0.4)^2))
  A <- (0.6*G1 + G2)*(vX + 1)
  
  B <- cos(12*sqrt(2*(x-0.4)^2 + (y-0.3)^2 + (z-0.2)^2)) * (sin(3*x)+0.2) * (cos(6*y) + 0.2) * (cos(8*z) + 0.2)
  
  C <- tanh(1+vX)*(1-z)
  
  return(R + (A + B + C + (vX + 2))/2)
}



###====================================
### Plotting True Function
###====================================


# Cube Texture Map
n_seq <- 201
x_seq <- seq(0,1, l=n_seq)
x_grid <- as.matrix(expand.grid(x_seq, x_seq))
xPl <- list(cbind(x_grid[,1], 0, x_grid[,2]), # front
            cbind(0, rev(x_grid[,2]), rev(x_grid[,1])), # left
            cbind(x_grid[,1], rev(x_grid[,2]), 0), # base
            cbind(1, rev(x_grid[,2]), x_grid[,1]), # right
            cbind(x_grid[,1], 1, rev(x_grid[,2])), # back
            cbind(x_grid[,1], x_grid[,2], 1)) # top

png(file="figures/T3D-f-cube-texture-map.png", type="cairo", width=3*1200, height=4*1200, res=400, bg="transparent")
par(mar=c(0,0,0,0), mfrow=c(4,3))
for (k in c(0, 1, 0, 2, 3, 4, 0, 5, 0, 0, 6, 0)){
  if(k==0){
    plot(0, 0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  } else {
    xP <- xPl[[k]]
    fP <- f(xP); fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)
    filled.contour3(x=x_seq, y=x_seq, z=fP_mat, color.palette=magma, levels=seq(-3.5, 3.5, 0.1), axes=F)
  }
}
dev.off()



# Plot Faces in Grid
n_seq <- 201
x_seq <- seq(0,1, l=n_seq)
x_grid <- as.matrix(expand.grid(x_seq, x_seq))


plot_f_slice <- function(xP, xlab, ylab, main=NULL) {
  fP <- f(xP)
  fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)
  filled.contour3(x=x_seq, y=x_seq, z=fP_mat, color.palette=magma, levels=seq(-3.5, 3.5, 0.1), xlab=xlab, ylab=ylab, main=main)
}



png(file="figures/T3D-f-grid.png", type="cairo", width=3*1200, height=3*1200, res=400, bg="transparent")
par(mar=c(4,4,2,2))
par(mfrow=c(3,3))
plot_f_slice(xP=cbind(0.5, x_grid[,1], x_grid[,2]), xlab="y", ylab="z", main="x=0.5")  # x=0.5
plot_f_slice(xP=cbind(x_grid[,1], x_grid[,2], z=1), xlab="x", ylab="y", main="z=1")  # z=1
plot_f_slice(xP=cbind(x_grid[,1], y=1, x_grid[,2]), xlab="x", ylab="z", main="y=1")  # y=1

plot_f_slice(xP=cbind(x_grid[,1], x_grid[,2], z=0), xlab="x", ylab="y", main="z=0")  # z=0
plot_f_slice(xP=cbind(x_grid[,1], 0.5, x_grid[,2]), xlab="x", ylab="z", main="y=0.5")  # y=0.5
plot_f_slice(xP=cbind(x=1, x_grid[,1], x_grid[,2]), xlab="y", ylab="z", main="x=1")  # x=1

plot_f_slice(xP=cbind(x_grid[,1], y=0, x_grid[,2]), xlab="x", ylab="z", main="y=0")  # y=0
plot_f_slice(xP=cbind(x=0, x_grid[,1], x_grid[,2]), xlab="y", ylab="z", main="x=0")  # x=0
plot_f_slice(xP=cbind(x_grid[,1], x_grid[,2], 0.5), xlab="x", ylab="y", main="z=0.5")  # z=0.5
dev.off()





###====================================
### Emulation - Design
###====================================

# From space-filling-designs
xD <- matrix(c(0, 31, 4,
               1, 62, 9,
               2, 13, 14,
               3, 44, 19,
               4, 75, 24,
               5, 26, 29,
               6, 57, 34,
               7, 8, 39,
               8, 39, 44,
               9, 70, 49,
               10, 21, 54,
               11, 52, 59,
               12, 3, 64,
               13, 34, 69,
               14, 65, 74,
               15, 16, 79,
               16, 47, 3,
               17, 78, 8,
               18, 29, 13,
               19, 60, 18,
               20, 11, 23,
               21, 42, 28,
               22, 73, 33,
               23, 24, 38,
               24, 55, 43,
               25, 6, 48,
               26, 37, 53,
               27, 68, 58,
               28, 19, 63,
               29, 50, 68,
               30, 1, 73,
               31, 32, 78,
               32, 63, 2,
               33, 14, 7,
               34, 45, 12,
               35, 76, 17,
               36, 27, 22,
               37, 58, 27,
               38, 9, 32,
               39, 40, 37,
               40, 71, 42,
               41, 22, 47,
               42, 53, 52,
               43, 4, 57,
               44, 35, 62,
               45, 66, 67,
               46, 17, 72,
               47, 48, 77,
               48, 79, 1,
               49, 30, 6,
               50, 61, 11,
               51, 12, 16,
               52, 43, 21,
               53, 74, 26,
               54, 25, 31,
               55, 56, 36,
               56, 7, 41,
               57, 38, 46,
               58, 69, 51,
               59, 20, 56,
               60, 51, 61,
               61, 2, 66,
               62, 33, 71,
               63, 64, 76,
               64, 15, 0,
               65, 46, 5,
               66, 77, 10,
               67, 28, 15,
               68, 59, 20,
               69, 10, 25,
               70, 41, 30,
               71, 72, 35,
               72, 23, 40,
               73, 54, 45,
               74, 5, 50,
               75, 36, 55,
               76, 67, 60,
               77, 18, 65,
               78, 49, 70,
               79, 0, 75), nrow = 80, ncol = 3, byrow = TRUE)/80 + 1/160

fD <- f(xD)


#plot_ly(x = xD[, 1], y = xD[, 2], z = xD[, 3], color = fD, type = "scatter3d", mode = "markers")



###====================================
### Embedding Surface Derivative
###====================================


s_deriv <- function(d) {
  (4*d * (0 <= d & d < 0.5)) + (-4*(d-1) * (0.5 <= d & d < 1))
}

v_deriv <- function(X){
  if (is.matrix(X)) { x <- X[,1]; y <- X[,2]; z <- X[,3] }
  else { x <- X[1]; y <- X[2]; z <- X[3] }
  
  #A <- 1 - s(3*sqrt((x-0.5)^2 + (0.5*z-0)^2 + 1e-6))
  dxA <- -s_deriv(3*sqrt((x-0.5)^2 + (0.5*z-0)^2)) * 3/(2*sqrt((x-0.5)^2 + (0.5*z-0)^2 + 1e-6)) * 2*(x-0.5)
  dyA <- rep(0, l=length(y))
  dzA <- -s_deriv(3*sqrt((x-0.5)^2 + (0.5*z-0)^2)) * 3/(2*sqrt((x-0.5)^2 + (0.5*z-0)^2 + 1e-6)) * 0.5*z
  
  B <- 1 - 2*(y > 0.5)
  
  dx <- dxA*B
  dy <- dyA*B
  dz <- dzA*B
  return(mapply(c, dx, dy, dz, SIMPLIFY=F))
}


# Design points fD, vD and derivatives
#vD <- v(xD)
#vD_deriv <- do.call(rbind, v_deriv(xD))

#plot_ly(x = xD[, 1], y = xD[, 2], z = xD[, 3], color = vD, type = "scatter3d", mode = "markers")
#plot_ly(x = xD[, 1], y = xD[, 2], z = xD[, 3], color = vD_deriv[,1], type = "scatter3d", mode = "markers")
#plot_ly(x = xD[, 1], y = xD[, 2], z = xD[, 3], color = vD_deriv[,2], type = "scatter3d", mode = "markers")
#plot_ly(x = xD[, 1], y = xD[, 2], z = xD[, 3], color = vD_deriv[,3], type = "scatter3d", mode = "markers")





###====================================
### TENSE
###====================================

#source("src/lib/TENSE.R")



CmatD <- k_NS(VX=V(xD), vdX=v_deriv(xD), k_S=k_S, theta=0.3, lambda2=1, sigma=1)







