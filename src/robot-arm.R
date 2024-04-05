library(viridisLite)

# Robot Arm Function `robot`
source("src/models/robot.R")

# Simple BL Emulator
source("src/lib/simple_BL_emulator.R")

# Maximin LHD (2D)
#source("src/lib/maximin_lhd_2d.R")

# Matrix wrapper for `robot` function
f <- function(X) {
  if (is.matrix(X)) {
    return(apply(X, 1, robot))
  } else {
    return(robot(X))
  }
}


var_cols <- function(n) hcl.colors(n, "YlOrRd", rev=TRUE)
surf_cols <- magma
diag_cols <- turbo


###=============================================
### Varying theta_2 and theta_4
###=============================================

ntheta <- 200
theta_seq <- seq(0, 2*pi, l=ntheta)

xP <- expand.grid(theta_seq, theta_seq)
n <- nrow(xP)
xP_all <- cbind(rep(0, n), xP[, 1], rep(pi/6, n), xP[, 2],
                rep(0.8, n), rep(0.4, n), rep(0.6, n), rep(1, n))

fP <- f(xP_all)
fP_mat <- matrix(fP, nrow=ntheta, ncol=ntheta)

tick_vals <- c(0, pi/2, pi, 3*pi/2, 2*pi)
tick_labs <- expression(0, pi/2, pi, 3*pi/2, 2*pi)

pdf(file="figures/2d-robot-arm-naive.pdf", width=7, height=6)
### True Function
filled.contour(x=theta_seq, y=theta_seq, z=fP_mat, color.palette=surf_cols, levels=seq(0, 3, 0.2),
               xlab=expression(theta[2]), ylab=expression(theta[4]),
               plot.axes={
                 axis(1, at=tick_vals, labels=tick_labs, las=1)
                 axis(2, at=tick_vals, labels=tick_labs, las=1)
                 contour(theta_seq, theta_seq, fP_mat, add=TRUE, level=seq(0, 3, 0.2),
                         lwd=0.4, drawlabels=FALSE)
               })


### Naive emulation without exploiting periodicity
#set.seed(6)
#xD <- maximin_lhd_2d(20, perturb=F)*2*pi
nD_seq <- 5
tD_seq <- seq(pi/nD_seq, 2*pi - pi/nD_seq, l=nD_seq)
xD <- expand.grid(tD_seq, tD_seq)
nD <- nrow(xD)
xD_all <- cbind(rep(0, nD), xD[, 1], rep(pi/6, nD), xD[, 2],
                rep(0.8, nD), rep(0.4, nD), rep(0.6, nD), rep(1, nD))

fD <- f(xD_all)


# Simple 2D BL emulator
BL_em <- simple_BL_emulator(xD=xD, D=fD, xP=xP, theta=2, sig=3, mu=1.5)
exp_mat <- matrix(BL_em[["ExpD_f(x)"]], nrow=ntheta, ncol=ntheta)
var_mat <- abs(matrix(BL_em[["VarD_f(x)"]], nrow=ntheta, ncol=ntheta))  # abs taken to correct for numerical instability

points_plot <- function(xD, col="black", cex=1.5, pch=20, ...) points(xD, col=col, cex=cex, pch=pch, ...)

# Emulator expectation plot
filled.contour(x=theta_seq, y=theta_seq, z=exp_mat, color.palette=surf_cols, levels=seq(0, 3, 0.2),
               xlab=expression(theta[2]), ylab=expression(theta[4]),
               plot.axes={
                 axis(1, at=tick_vals, labels=tick_labs, las=1)
                 axis(2, at=tick_vals, labels=tick_labs, las=1)
                 contour(theta_seq, theta_seq, exp_mat, add=TRUE, level=seq(0, 3, 0.2),
                         lwd=0.4, drawlabels=FALSE)
                 points_plot(xD)
               })

filled.contour(x=theta_seq, y=theta_seq, z=exp_mat, color.palette=surf_cols, levels=seq(0, 3, 0.2),
               xlab=expression(theta[2]), ylab=expression(theta[4]),
               plot.axes={
                 axis(1, at=tick_vals, labels=tick_labs, las=1)
                 axis(2, at=tick_vals, labels=tick_labs, las=1)
                 contour(theta_seq, theta_seq, exp_mat, add=TRUE, level=seq(0, 3, 0.2),
                         lwd=0.4, drawlabels=FALSE)
               })

# Emulator sd plot
filled.contour(x=theta_seq, y=theta_seq, z=sqrt(var_mat), color.palette=var_cols, levels=seq(0, 1, 0.1),
               xlab=expression(theta[2]), ylab=expression(theta[4]),
               plot.axes={
                 axis(1, at=tick_vals, labels=tick_labs, las=1)
                 axis(2, at=tick_vals, labels=tick_labs, las=1)
                 points_plot(xD)
               })

# Emulator diagnostics plot
diag_mat <- (exp_mat - fP_mat) / sqrt(var_mat)
filled.contour(x=theta_seq, y=theta_seq, z=diag_mat, color.palette=diag_cols, levels=seq(-3,3,0.25),
               xlab=expression(theta[2]), ylab=expression(theta[4]),
               plot.axes={axis(1);axis(2)
                 contour(theta_seq, theta_seq, diag_mat, add=TRUE, level=seq(-3, 3, 0.25),
                         lwd=0.4, drawlabels=FALSE)
                 points_plot(xD)
               })
dev.off()



### Demonstrating periodicity

ntheta <- 400
theta_seq <- seq(-pi, 3*pi, l=ntheta)

grid <- expand.grid(theta_seq, theta_seq)
n <- nrow(grid)
xP <- cbind(rep(0, n), grid[, 1], rep(pi/6, n), grid[, 2],
            rep(0.8, n), rep(0.4, n), rep(0.6, n), rep(1, n))

fP <- f(xP)
fP_mat <- matrix(fP, nrow=ntheta, ncol=ntheta)

tick_vals <- c(0, 2*pi)
tick_labs <- expression(0, 2*pi)


### True function with extended ranges
pdf(file="figures/2d-robot-arm-periodicity.pdf", width=7, height=6)
filled.contour(x=theta_seq, y=theta_seq, z=fP_mat, color.palette=surf_cols, levels=seq(0, 2.8, 0.1),
               xlab=expression(theta[2]), ylab=expression(theta[4]),
               plot.axes={
                 axis(1, at=tick_vals, labels=tick_labs, las=1)
                 axis(2, at=tick_vals, labels=tick_labs, las=1)
                 lines(c(0, 0),        c(-2*pi, 4*pi), lwd=2)
                 lines(c(2*pi, 2*pi),  c(-2*pi, 4*pi), lwd=2)
                 lines(c(-2*pi, 4*pi), c(0, 0),        lwd=2)
                 lines(c(-2*pi, 4*pi), c(2*pi, 2*pi),  lwd=2)
                 box(lwd=2)
                 lines(c(0, 2*pi),    c(0, 0),       lwd=4)
                 lines(c(0, 2*pi),    c(2*pi, 2*pi), lwd=4)
                 lines(c(0, 0),       c(0, 2*pi),    lwd=4)
                 lines(c(2*pi, 2*pi), c(0, 2*pi),    lwd=4)
               })
dev.off()

