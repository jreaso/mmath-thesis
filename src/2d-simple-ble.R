library(viridisLite)

# Simple BL Emulator
source("src/lib/simple_BL_emulator.R")


### Plotting setup
axis_labels <- c("", "")  #c(expression(x), expression(y))
points_plot <- function(xD, col="black", cex=1.5, pch=20, ...) points(xD, col=col, cex=cex, pch=pch, ...)

var_cols <- function(n) hcl.colors(n, "YlOrRd", rev=TRUE)
surf_cols <- magma
diag_cols <- turbo


### True Function
f <- function(X){
  A <- exp(-10*((X[,1]-0.4)^2 + (X[,2]-0.7)^2))
  B <- exp(-20*((X[,1]-0.9)^2 + (X[,2]-0.5)^2))
  C <- cos(6*(X[,1]^2 - X[,2]^2))
  return(A + 5/3 * B + 1/3 * C - 0.5)
}


### Plotting
n_seq <- 201
x_seq <- seq(0,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))

# Evaluate true function
fP <- f(xP)
fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)  # resize to matrix for plotting

# Plot true function
pdf(file="figures/2d-simple-ble.pdf", width=7, height=6)
par(mar=c(4,4,2,4))
filled.contour(x=x_seq, y=x_seq, z=fP_mat, color.palette=surf_cols, levels=seq(-1, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, fP_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
               })


nD_seq <- 4
xD_seq <- seq(0.05, 0.95, len=nD_seq)

xD <- as.matrix(expand.grid(xD_seq, xD_seq))
fD <- f(xD)

# Simple 2D BL emulator
BL_em <- simple_BL_emulator(xD=xD, D=fD, xP=xP, theta=0.35, sig=1, mu=0)
exp_mat <- matrix(BL_em[["ExpD_f(x)"]], nrow=n_seq, ncol=n_seq)
var_mat <- matrix(BL_em[["VarD_f(x)"]], nrow=n_seq, ncol=n_seq)


# Emulator expectation plot
filled.contour(x=x_seq, y=x_seq, z=exp_mat, color.palette=surf_cols, level=seq(-1, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, exp_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
               })

# Emulator variance plot
filled.contour(x=x_seq, y=x_seq, z=var_mat, color.palette=var_cols, levels=seq(0,0.1,0.01),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 points_plot(xD)
               })

# Emulator sd plot
sd_mat <- sqrt(abs(var_mat)) # abs taken to avoid NaNs from error with decimal numbers close to 0
filled.contour(x=x_seq, y=x_seq, z=sd_mat, color.palette=var_cols, levels=seq(0,0.3,0.025),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 points_plot(xD)
               })

# Emulator diagnostics plot
diag_mat <- (exp_mat - fP_mat) / sd_mat
filled.contour(x=x_seq, y=x_seq, z=diag_mat, color.palette=diag_cols, levels=seq(-3,3,0.25),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, diag_mat, add=TRUE, level=seq(-3, 3, 0.25), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD)
               })
dev.off()

