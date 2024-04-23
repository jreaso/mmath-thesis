library(viridisLite)

source("src/lib/simple_BL_emulator.R")
source("src/lib/maximin_lhd_2d.R")
source("src/lib/sequential_minimax.R")

surf_cols <- magma
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo

axis_labels <- c(expression(x), expression(y))
points_plot <- function(xD, col="black", bg="black", cex=1.5, pch=21, ...) points(xD, col=col, cex=cex, bg=bg, pch=pch, ...)


f <- function(x) tanh((sin(3*pi*x[,1])+cos(3*pi*x[,1]*x[,2]))*exp(-3*((x[,1]-0.5)^2+(x[,2]-0.5)^2)))

n_seq <- 101
x_seq <- seq(0, 1, len=n_seq)
xP <- as.matrix(expand.grid("x1"=x_seq, "x2"=x_seq))


fP <- f(xP)
fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)


pdf(file="figures/2d-history-match.pdf", width=7, height=6)
par(mar=c(4,4,2,4))

# True Function
filled.contour(x=x_seq, y=x_seq, z=fP_mat, color.palette=surf_cols, levels=seq(-1.6, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, fP_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
               })


### History Matching
z <- 0.85
sigma_e <- 0.02
sigma_epsilon <- 0.04


imp_cols <- function(n) turbo(n, begin=0.15, end=1)
imp_levs <- c(0, seq(1, 2.75, 0.25), seq(3, 15, 2), 20, 30, 50)

Imp_true_mat <- sqrt( (fP_mat - z)^2 / (sigma_e^2 + sigma_epsilon^2) ) 

# True Implausibility
filled.contour(x=x_seq, y=x_seq, z=Imp_true_mat, color.palette=imp_cols, levels=imp_levs,
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 #contour(x_seq, x_seq, Imp_true_mat, add=TRUE, level=imp_levs, lwd=0.4, drawlabels=FALSE)
                 contour(x_seq, x_seq, Imp_true_mat, add=TRUE, level=c(3), lwd=2, drawlabels=T)
               })






######## xD points to choose from
#20x20 grid of points to choose from for design
xD_seq <- seq(1/40, 1-1/40, len=20)
xD_base <- round(as.matrix(expand.grid("x1"=xD_seq, "x2"=xD_seq)), 3)




###### WAVE 1


set.seed(1)
# LHD
xD_w1 <- maximin_lhd_2d(nl=20, perturb=F)
xD_w1 <- round(xD_w1, 3)
fD1 <- f(xD_w1)


BL_em1 <- simple_BL_emulator(xD=xD_w1, D=fD1, xP=xP, theta=0.28, sig=1, mu=0)
exp1_mat <- matrix(BL_em1[["ExpD_f(x)"]], nrow=n_seq, ncol=n_seq)
var1_mat <- matrix(BL_em1[["VarD_f(x)"]], nrow=n_seq, ncol=n_seq)


# wave 1 expectation
filled.contour(x=x_seq, y=x_seq, z=exp1_mat, color.palette=surf_cols, level=seq(-1.6, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, exp1_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD_w1, bg="white")
               })

# wave 1 standard deviation
sd1_mat <- sqrt(abs(var1_mat))
filled.contour(x=x_seq, y=x_seq, z=sd1_mat, color.palette=var_cols, levels=seq(0,1,0.05),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 points_plot(xD_w1, bg="white")
               })


# Emulator diagnostics plot
diag1_mat <- (exp1_mat - fP_mat) / sd1_mat
filled.contour(x=x_seq, y=x_seq, z=diag1_mat, color.palette=diag_cols, levels=seq(-3,3,0.25),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, diag1_mat, add=TRUE, level=seq(-3, 3, 0.25), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD_w1, bg="white")
               })


### Implausibility
Imp1_mat <- sqrt((exp1_mat - z)^2 / (var1_mat + sigma_e^2 + sigma_epsilon^2))
filled.contour(x=x_seq, y=x_seq, z=Imp1_mat, color.palette=imp_cols, levels=imp_levs,
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 #contour(x_seq, x_seq, Imp1_mat, add=TRUE, level=imp_levs, lwd=0.4, drawlabels=FALSE)
                 contour(x_seq, x_seq, Imp1_mat, add=TRUE, level=c(3), lwd=2, drawlabels=T)
                 points_plot(xD_w1, bg="white")
                 points_plot(xD_w1, col="black", pch=20, cex=1)
               })


BL_em_base <- simple_BL_emulator(xD=xD_w1, D=fD1, xP=xD_base, theta=0.28, sig=1, mu=0)
exp_base <- BL_em_base[["ExpD_f(x)"]]
var_base <- BL_em_base[["VarD_f(x)"]]

imp_base <- sqrt((exp_base - z)^2 / (var_base + sigma_e^2 + sigma_epsilon^2))

xD_base_w1 <- xD_base[imp_base <= 3, ]




###### WAVE 2


set.seed(42)
xD_w1_inds <- as.vector(na.omit(sapply(1:nrow(xD_w1), function(i) {
  which(xD_base_w1[, 1] == xD_w1[i, 1] & xD_base_w1[, 2] == xD_w1[i, 2])[1]
})))
xD_w2_inds <- sequential_minimax(X=xD_base_w1, k=12+length(xD_w1_inds), p=3, niter=100000, XDi=xD_w1_inds)
xD_w2 <- xD_base_w1[setdiff(xD_w2_inds, xD_w1_inds), ]

xD_w2_full <- rbind(xD_w1, xD_w2)

fD2 <- f(xD_w2_full)


BL_em2 <- simple_BL_emulator(xD=xD_w2_full, D=fD2, xP=xP, theta=0.28, sig=1, mu=0)
exp2_mat <- matrix(BL_em2[["ExpD_f(x)"]], nrow=n_seq, ncol=n_seq)
var2_mat <- matrix(BL_em2[["VarD_f(x)"]], nrow=n_seq, ncol=n_seq)


# wave 1 expectation
filled.contour(x=x_seq, y=x_seq, z=exp2_mat, color.palette=surf_cols, level=seq(-1.6, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, exp2_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD_w2_full)
                 points_plot(xD_w2, bg="white")
               })

# wave 1 standard deviation
sd2_mat <- sqrt(abs(var2_mat))
filled.contour(x=x_seq, y=x_seq, z=sd2_mat, color.palette=var_cols, levels=seq(0,1,0.05),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 points_plot(xD_w2_full)
                 points_plot(xD_w2, bg="white")
               })


# Emulator diagnostics plot
diag2_mat <- (exp2_mat - fP_mat) / sd2_mat
filled.contour(x=x_seq, y=x_seq, z=diag2_mat, color.palette=diag_cols, levels=seq(-3,3,0.25),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, diag2_mat, add=TRUE, level=seq(-3, 3, 0.25), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD_w2_full)
                 points_plot(xD_w2, bg="white")
               })


### Implausibility
Imp2_mat <- sqrt((exp2_mat - z)^2 / (var2_mat + sigma_e^2 + sigma_epsilon^2))
filled.contour(x=x_seq, y=x_seq, z=Imp2_mat, color.palette=imp_cols, levels=imp_levs,
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 #contour(x_seq, x_seq, Imp2_mat, add=TRUE, level=imp_levs, lwd=0.4, drawlabels=FALSE)
                 contour(x_seq, x_seq, Imp2_mat, add=TRUE, level=c(3), lwd=2, drawlabels=T)
                 points_plot(xD_w2_full, bg="white")
                 points_plot(xD_w2, col="black", pch=20, cex=1)
               })


BL_em_base <- simple_BL_emulator(xD=xD_w2_full, D=fD2, xP=xD_base, theta=0.28, sig=1, mu=0)
exp_base <- BL_em_base[["ExpD_f(x)"]]
var_base <- BL_em_base[["VarD_f(x)"]]

imp_base <- sqrt((exp_base - z)^2 / (var_base + sigma_e^2 + sigma_epsilon^2))

xD_base_w2 <- xD_base[imp_base <= 3, ]




###### WAVE 3


set.seed(42)
xD_w2_inds <- as.vector(na.omit(sapply(1:nrow(xD_w2_full), function(i) {
  which(xD_base_w2[, 1] == xD_w2_full[i, 1] & xD_base_w2[, 2] == xD_w2_full[i, 2])[1]
})))
xD_w3_inds <- sequential_minimax(X=xD_base_w2, k=20+length(xD_w2_inds), p=3, niter=100000, XDi=xD_w2_inds)
xD_w3 <- xD_base_w2[setdiff(xD_w3_inds, xD_w2_inds), ]

xD_w3_full <- rbind(xD_w3, xD_w2_full)

fD3 <- f(xD_w3_full)


BL_em3 <- simple_BL_emulator(xD=xD_w3_full, D=fD3, xP=xP, theta=0.28, sig=1, mu=0)
exp3_mat <- matrix(BL_em3[["ExpD_f(x)"]], nrow=n_seq, ncol=n_seq)
var3_mat <- matrix(BL_em3[["VarD_f(x)"]], nrow=n_seq, ncol=n_seq)


# wave 1 expectation
filled.contour(x=x_seq, y=x_seq, z=exp3_mat, color.palette=surf_cols, level=seq(-1.6, 1, 0.1),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, exp3_mat, add=TRUE, level=seq(-1, 1, 0.1), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD_w3_full)
                 points_plot(xD_w3, bg="white")
               })

# wave 1 standard deviation
sd3_mat <- sqrt(abs(var3_mat))
filled.contour(x=x_seq, y=x_seq, z=sd3_mat, color.palette=var_cols, levels=seq(0,1,0.05),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 points_plot(xD_w3_full)
                 points_plot(xD_w3, bg="white")
               })


# Emulator diagnostics plot
diag3_mat <- (exp3_mat - fP_mat) / sd3_mat
filled.contour(x=x_seq, y=x_seq, z=diag3_mat, color.palette=diag_cols, levels=seq(-3,3,0.25),
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 contour(x_seq, x_seq, diag3_mat, add=TRUE, level=seq(-3, 3, 0.25), lwd=0.4, drawlabels=FALSE)
                 points_plot(xD_w3_full)
                 points_plot(xD_w3, bg="white")
               })


### Implausibility
Imp3_mat <- sqrt((exp3_mat - z)^2 / (var3_mat + sigma_e^2 + sigma_epsilon^2))
filled.contour(x=x_seq, y=x_seq, z=Imp3_mat, color.palette=imp_cols, levels=imp_levs,
               xlab=axis_labels[1], ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 #contour(x_seq, x_seq, Imp3_mat, add=TRUE, level=imp_levs, lwd=0.4, drawlabels=FALSE)
                 contour(x_seq, x_seq, Imp3_mat, add=TRUE, level=c(3), lwd=2, drawlabels=T)
                 points_plot(xD_w3_full, bg="white")
                 points_plot(xD_w3, col="black", pch=20, cex=1)
               })


BL_em_base <- simple_BL_emulator(xD=xD_w3_full, D=fD3, xP=xD_base, theta=0.28, sig=1, mu=0)
exp_base <- BL_em_base[["ExpD_f(x)"]]
var_base <- BL_em_base[["VarD_f(x)"]]

imp_base <- sqrt((exp_base - z)^2 / (var_base + sigma_e^2 + sigma_epsilon^2))

xD_base_w3 <- xD_base[imp_base <= 3, ]


dev.off()

