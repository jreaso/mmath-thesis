library(viridisLite)
library(plotly)
library(plot3D)

source("src/lib/TENSE.R")
source("src/lib/advanced_NS_emulator.R")
source("src/lib/filled.contour3.R")
source("src/lib/simple_NS_emulator2.R")

var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo
surf_cols <- magma #plasma


###=========================================
### Embedding Surface
###=========================================

# Sigmoid like function (derivative has a compact domain on [0, 1])
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

### Derivatives
s_deriv <- function(d) {
  (4*d * (0 <= d & d < 0.5)) + (-4*(d-1) * (0.5 <= d & d < 1))
}

v_deriv <- function(X){
  if (is.matrix(X)) { x <- X[,1]; y <- X[,2]; z <- X[,3] }
  else { x <- X[1]; y <- X[2]; z <- X[3] }
  
  #A <- 1 - s(3*sqrt((x-0.5)^2 + (0.5*z-0)^2 + 1e-6))
  dxA <- -s_deriv(3*sqrt((x-0.5)^2 + (0.5*z-0)^2)) * 3/(2*sqrt((x-0.5)^2 + (0.5*z-0)^2 + 1e-6)) * 2*(x-0.5)
  dzA <- -s_deriv(3*sqrt((x-0.5)^2 + (0.5*z-0)^2)) * 3/(2*sqrt((x-0.5)^2 + (0.5*z-0)^2 + 1e-6)) * 0.5*z
  
  B <- 1 - 2*(y > 0.5)
  
  dx <- dxA*B
  dz <- dzA*B
  return(cbind(dx, 0, dz))
}


### Plotting Embedding Surface

# Compact domain sigmoid plot
pdf(file="figures/TENSE-3D/1d-sigmoidal.pdf", width=6, height=5)
par(mar = c(4, 4, 2, 2))
plot(seq(-0.5, 1.5, l=100), s(seq(-0.5, 1.5, l=100)), type="l", xlab="d", ylab="s(d)", xlim=c(-0.1,1.1))
dev.off()


n_seq <- 101
x_seq <- seq(0,1, l=n_seq);
x_grid <- as.matrix(expand.grid(x_seq, x_seq))

# Slice (a)
xP <- cbind(x_grid[,1], 0.8, x_grid[,2])
vP <- v(xP); vP_mat <- matrix(vP, nrow=n_seq, ncol=n_seq)
png(file="figures/TENSE-3d/v-a.png", type="cairo", width=2400, height=2400, res=600, bg="transparent")
par(mar=c(4,4,2,2))
filled.contour3(x=x_seq, y=x_seq, z=vP_mat, color.palette=magma, levels=seq(-1, 1, 0.05), xlab="x", ylab="z")
dev.off()

# Slice (b)
xP <- cbind(x_grid[,1], 0.2, x_grid[,2])
vP <- v(xP); vP_mat <- matrix(vP, nrow=n_seq, ncol=n_seq)
png(file="figures/TENSE-3d/v-b.png", type="cairo", width=2400, height=2400, res=600, bg="transparent")
par(mar=c(4,4,2,2))
filled.contour3(x=x_seq, y=x_seq, z=vP_mat, color.palette=magma, levels=seq(-1, 1, 0.05), xlab="x", ylab="z")
dev.off()

# Slice (c)
xP <- cbind(0.5, x_grid)
vP <- v(xP); vP_mat <- matrix(vP, nrow=n_seq, ncol=n_seq)
png(file="figures/TENSE-3d/v-c.png", type="cairo", width=2400, height=2400, res=600, bg="transparent")
par(mar=c(4,4,2,2))
filled.contour3(x=x_seq, y=x_seq, z=vP_mat, color.palette=magma, levels=seq(-1, 1, 0.05), xlab="y", ylab="z",
                plot.axes={axis(1);axis(2)
                  lines(c(0.5, 0.5), c(0, 2/3), type="l", lwd=3)
                })
dev.off()

# Legend
png(file="figures/TENSE-3d/v-legend.png", type="cairo", width=900, height=2400, res=600, bg="transparent")
par(mar=c(4,0,2,2))
colkey(col=magma(40), side=4, clim=c(-1, 1), add=F, width=4)
dev.off()


###=========================================
### True Function
###=========================================

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


### Plotting True Function

n_seq <- 101
x_seq <- seq(0,1, l=n_seq);
x_grid <- as.matrix(expand.grid(x_seq, x_seq))

names <- c("face-front", "face-back", "face-base", "face-top", "face-left", "face-right",
           "slice-xz1", "slice-xz2", "slice-yz")

xP_list <- list(cbind(x_grid[,1], 0, x_grid[,2]), cbind(x_grid[,1], 1, x_grid[,2]),
                cbind(x_grid, 0), cbind(x_grid, 1), cbind(0, x_grid), cbind(1, x_grid),
                cbind(x_grid[,1], 0.2, x_grid[,2]), cbind(x_grid[,1], 0.8, x_grid[,2]),
                cbind(0.5, x_grid))

mains <- c("y=0", "y=1", "z=0", "z=1", "x=0", "x=1", "y=0.2", "y=0.8", "x=0.5")
xlabs <- c("x", "x", "x", "x", "y", "y", "x", "x", "y")
ylabs <- c("z", "z", "y", "y", "z", "z", "z", "z", "z")

disc_plots <- list(NULL, NULL,
                   function() lines(c(1/6, 5/6), c(0.5, 0.5), type="l", lwd=2), # face-base
                   NULL, NULL, NULL, NULL, NULL,
                   function() lines(c(0.5, 0.5), c(0, 0.6), type="l", lwd=2)) # slice-yz


# True Function Plots
for (i in 1:9){
  xP <- xP_list[[i]]
  fP <- f(xP)
  fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)
  png(file=paste0("figures/TENSE-3D/f-", names[i], ".png"), type="cairo", width=2400, height=2400, res=600, bg="transparent")
  par(mar=c(4,4,2,2), mfrow=c(1,1))
  filled.contour3(x=x_seq, y=x_seq, z=fP_mat, color.palette=surf_cols, levels=seq(-3.5, 3.5, 0.1),
                  xlab=xlabs[i], ylab=ylabs[i], #main=mains[i],
                  plot.axes={axis(1);axis(2)
                    if (!is.null(disc_plots[[i]])){
                      disc_plots[[i]]()
                    }
                  })
  dev.off()
}



# Legend
png(file="figures/TENSE-3d/surf-legend.png", type="cairo", width=900, height=2400, res=600, bg="transparent")
par(mar=c(4,0,2,2))
colkey(col=surf_cols(70), side=4, clim=c(-3.5, 3.5), add=F, width=4)
dev.off()




###=========================================
### TENSE Emulation
###=========================================
#source("src/lib/TENSE.R")
#source("src/lib/advanced_NS_emulator.R")


### Design
load("data/designs/Mm3_160.RData")
xD <- Mm3_160/160 + 1/320
colnames(xD) <- c("x", "y", "z")

fD <- f(xD)
plot_ly(x = xD[, 1], y = xD[, 2], z = xD[, 3], color = fD, type = "scatter3d", mode = "markers")

pt_cexs <- list(1-abs(xD[,2]), 1-abs(xD[,2]-1),
                1-abs(xD[,3]), 1-abs(xD[,3]-1), 1-abs(xD[,1]), 1-abs(xD[,1]-1),
                1-abs(xD[,2]-0.2), 1-abs(xD[,2]-0.8),
                1-abs(xD[,1]-0.5))


# Optionally use smaller grid for emulation
# WARNING: Emulation takes a while with a large grid size 
#n_seq <- 21
#x_grid <- as.matrix(expand.grid(x_seq, x_seq))
#x_seq <- seq(0,1, l=n_seq);
#xP_list <- list(cbind(x_grid[,1], 0, x_grid[,2]), cbind(x_grid[,1], 1, x_grid[,2]),
#                cbind(x_grid, 0), cbind(x_grid, 1), cbind(0, x_grid), cbind(1, x_grid),
#                cbind(x_grid[,1], 0.2, x_grid[,2]), cbind(x_grid[,1], 0.8, x_grid[,2]),
#                cbind(0.5, x_grid))


em_out_list <- list()


### Advanced Emulator
VD <- V(xD)
fD <- f(xD)
vdxD <- v_deriv(xD) 

# Basis Functions
# 1, x, y, z, x^2, y^2, z^2, xz
g <- list(function(X) 1, function(X) X[1], function(X) X[2], function(X) X[3], function(X) X[1]^2, function(X) X[2]^2, function(X) X[3]^2, function(X) X[1]*X[3])
mu_beta <- c(-2.5, -4, 0, -1.5, 3, 2, 0, 5)
Sigma_beta <- 10*diag(length(g))

for (i in 1:9){
  xP <- xP_list[[i]]
  VP <- V(xP)
  vdxP <- v_deriv(xP)
  em_out <- advanced_NS_emulator(VxD=VD, D=fD, VxP=VP, vdxD=vdxD, vdxP=vdxP, k_S=k_S, mu_beta=mu_beta, g=g, Sigma_beta=Sigma_beta, theta=0.22, lambda2=0.06, sigma=1, just_var=T)
  # Simple TENSE Emulator
  #simple_NS_emulator2(VxD=VD, D=fD, VxP=VP, vdxD=vdxD, vdxP=vdxP, k_S=k_S, mu=0, theta=0.22, lambda2=0.06, sigma=1, just_var=T)
  em_out_list[[i]] <- em_out
  print(paste0("Completed ",i,"/9"))
}



### Regression Coefficients
reg_coefs <- as.vector(em_out_list[[1]][["ED_beta"]])
quants <- c("1", "x", "y", "z", "x^2", "y^2", "z^2", "xz")
for (j in 1:8) {
  print(paste0(quants[j], ": ", round(reg_coefs[j], 3)))
}


### Emulator Plots
for (i in 1:9){
  em_out <- em_out_list[[i]]
  
  # Expectation
  exp <- em_out[["ED_fP"]]
  exp_mat <- matrix(exp, nrow=n_seq, ncol=n_seq)
  
  png(file=paste0("figures/TENSE-3D/exp-", names[i], ".png"), type="cairo", width=2400, height=2400, res=600, bg="transparent")
  par(mar=c(4,4,2,2), mfrow=c(1,1))
  filled.contour3(x=x_seq, y=x_seq, z=exp_mat, color.palette=surf_cols, levels=seq(-3.5, 3.5, 0.1),
                  xlab=xlabs[i], ylab=ylabs[i], #main=mains[i],
                  plot.axes={axis(1);axis(2)
                    if (!is.null(disc_plots[[i]])){
                      disc_plots[[i]]()
                    }
                  })
  dev.off()
  
  
  # Standard Deviation
  var <- em_out[["VarD_fP"]] + 1e-10 # for stability
  var_mat <- matrix(var, nrow=n_seq, ncol=n_seq)
  sd_mat <- sqrt(var_mat)
  
  png(file=paste0("figures/TENSE-3D/sd-", names[i], ".png"), type="cairo", width=2400, height=2400, res=600, bg="transparent")
  filled.contour3(x=x_seq, y=x_seq, z=sd_mat, color.palette=var_cols, levels=seq(0, 1.5, 0.05),
                  xlab=xlabs[i], ylab=ylabs[i], #main=mains[i],
                  plot.axes={axis(1);axis(2)
                    if (!is.null(disc_plots[[i]])){
                      disc_plots[[i]]()
                    }
                    #points(x=xD[, xlabs[i]], y=xD[, ylabs[i]], pch=21, bg=1, cex=pt_cexs[[i]])
                  })
  dev.off()
  
  
  # Diagnostics
  xP <- xP_list[[i]]
  fP <- f(xP)
  fP_mat <- matrix(fP, nrow=n_seq, ncol=n_seq)
  diag_mat <- (exp_mat - fP_mat)/sd_mat
  
  png(file=paste0("figures/TENSE-3D/diag-", names[i], ".png"), type="cairo", width=2400, height=2400, res=600, bg="transparent")
  par(mar=c(4,4,2,2), mfrow=c(1,1))
  filled.contour3(x=x_seq, y=x_seq, z=diag_mat, color.palette=diag_cols, levels=seq(-3, 3, 0.25),
                  xlab=xlabs[i], ylab=ylabs[i], #main=mains[i],
                  plot.axes={axis(1);axis(2)
                    contour(x_seq, x_seq, diag_mat, add=TRUE, level=seq(-3, 3, 0.25), lwd=0.4, drawlabels=FALSE)
                    if (!is.null(disc_plots[[i]])){
                      disc_plots[[i]]()
                    }
                  })
  dev.off()
}


# Legend
png(file="figures/TENSE-3d/var-legend.png", type="cairo", width=900, height=2400, res=600, bg="transparent")
par(mar=c(4,0,2,2))
colkey(col=var_cols(20), side=4, clim=c(0, 1.5), add=F, width=4)
dev.off()

# Legend
png(file="figures/TENSE-3d/diag-legend.png", type="cairo", width=900, height=2400, res=600, bg="transparent")
par(mar=c(4,0,2,2))
colkey(col=diag_cols(24), side=4, clim=c(-3, 3), add=F, width=4)
dev.off()








# Optionally can optimise choice of theta and lambda2
#n_seq <- 21
#x_seq <- seq(0,1, l=n_seq);
#x2_grid <- as.matrix(expand.grid(x_seq, x_seq))
#xP <- cbind(0.5, x2_grid[,1], x2_grid[,2])

#VD <- V(xD); VP <- V(xP)
#fD <- f(xD)
#vdxD <- v_deriv(xD); vdxP <- v_deriv(xP)
#fP <- f(xP)

#wrapper <- function(params) {
#  theta <- params[1]
#  lambda2 <- abs(params[2])
#  
#  exp <- simple_NS_emulator2(VxD=VD, D=fD, VxP=VP, vdxD=vdxD, vdxP=vdxP, k_S=k_S, mu=0, theta=theta, lambda2=lambda2, sigma=1, just_var=T)[["ED_fP"]]
#  return(mean((exp-fP)^2)) # MSE
#}

#result <- optim(par=c(0.3, 0.1), fn=wrapper, method="L-BFGS-B") #"L-BFGS-B" or "Nelder-Mead" (simplex method)
#result$par




