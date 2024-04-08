library(LHD)
library(lhs)
library(viridisLite)

# Simple BL Emulator
source("src/lib/simple_BL_emulator.R")

# 2D maximin LHD
source("src/lib/maximin_lhd_2d.R")


### Plotting setup
axis_labels <- c("", "")  #c(expression(x), expression(y))
points_plot <- function(xD, col="black", bg="springgreen2", cex=1.5, pch=21, ...) points(xD, col=col, bg=bg, cex=cex, pch=pch, ...)

var_cols <- function(n) hcl.colors(n, "YlOrRd", rev=TRUE)
surf_cols <- magma
diag_cols <- turbo


# Test function 1
f <- function(X) {
  if (is.matrix(X)) {
    x <- X[,1]; y <- X[,2]
  } else {
    x <- X[1]; y <- X[2]
  }
  
  return(1.3*(sin(4*pi*x)^2 - 0.4 + (x^2 * y)/10 - (y-0.5)^2))
}


n_seq <- 501
x_seq <- seq(0,1, len=n_seq)
xP <- as.matrix(expand.grid(x_seq, x_seq))

fP <- f(xP)
fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)


nD_seq <- 4
xD_seq <- seq(1/(2*nD_seq), 1-1/(2*nD_seq), len=nD_seq)

xD <- as.matrix(expand.grid(xD_seq, xD_seq))
fD <- f(xD)


# Plot true function
pdf(file="figures/2d-grid-design-flaws.pdf", width=7, height=6)
par(mar=c(4,4,2,4))
filled.contour(x=x_seq, y=x_seq, z=fP_mat, color.palette=surf_cols, levels=seq(-1, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, fP_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD)
               })



# Simple 2D BL emulator
BL_em <- simple_BL_emulator(xD=xD, D=fD, xP=xP, theta=0.35, sig=1, mu=0)
exp_mat <- matrix(BL_em[["ExpD_f(x)"]], nrow=n_seq, ncol=n_seq)
var_mat <- abs(matrix(BL_em[["VarD_f(x)"]], nrow=n_seq, ncol=n_seq))
# abs taken to avoid NaNs from error with decimal numbers close to 0


# Emulator expectation plot
filled.contour(x=x_seq, y=x_seq, z=exp_mat, color.palette=surf_cols, level=seq(-1, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, exp_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD)
               })

# Emulator var plot
filled.contour(x=x_seq, y=x_seq, z=var_mat, color.palette=var_cols, levels=seq(0,0.25,0.025),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 points_plot(xD)
               })

# Emulator diagnostics plot
sd_mat <- sqrt(var_mat)
diag_mat <- (exp_mat - fP_mat) / sd_mat
filled.contour(x=x_seq, y=x_seq, z=diag_mat, color.palette=diag_cols, levels=seq(-3,3,0.25),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, diag_mat, add=TRUE, level=seq(-3, 3, 0.25), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD)
               })
dev.off()


