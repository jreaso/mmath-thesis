library(viridisLite)

source("src/lib/simple_BL_emulator.R")

f <- function(x) {
  sin(2*pi*x + 1) + 0.6*cos(5*pi*x) 
}

xP <- seq(-0.05, 1.05, len=1001)
fxP <- lapply(xP, f)

xD <- seq(0.05,0.95, len=6)
D <- f(xD)

# Simple Emulator
em_out <- simple_BL_emulator(xD=as.matrix(xD, ncol=1), D=D, xP=as.matrix(xP, ncol=1),
                             theta=.25, sig=1, nugget=0, mu=0)
em_exp <- as.vector(em_out[["ExpD_f(x)"]])
em_var <- as.vector(em_out[["VarD_f(x)"]])

int_col <- rgb(.9,.9,.9)

pdf(file="figures/1d-implausibility.pdf", width=6, height=4.5)
plot(1, 1, xlim=c(0,1), ylim=c(-3,3), type="n", xlab="", ylab="")
upper <- em_exp + 3*sqrt(em_var); lower <- em_exp - 3*sqrt(em_var)
polygon(c(xP, rev(xP)), c(upper, rev(lower)), col=int_col, border=NA)
#lines(xP, fxP, col='red', lwd=2.5)
lines(xP, em_exp, lwd=2, col=1)#, lty=2.5)
points(xD, D, pch=21, col=1, bg=1, cex=1)
box()
title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)


z <- -1.2
sigma_e <- 0.08
sigma_eps <- 0.05
v <- (sigma_e^2+ sigma_eps^2)
z_error <- c(z-3*sqrt(v), z+3*sqrt(v))

imp <- sqrt((em_exp - z)^2 / (em_var + v))


# Implausibility Plot
plot(xP, imp, type="l", lwd=2, xlab="x", ylab="I(x)", xlim=c(0,1), ylim=c(-0.9,30))
abline(h=3, lwd=2, col="red")

imp_breaks <- c(seq(0,2.8,0.2),seq(3,9.5,0.5),seq(10,80,10))
imp_cut_index <- cut(imp, imp_breaks, labels=FALSE)
cols <- turbo(length(imp_breaks)-1,begin=0.15,end=1)
rect_wid <- 1.1/(length(xP)-1)
for(i in 1:length(xP)) rect(xP[i]-rect_wid/2,-2.2,xP[i]+rect_wid/2, 0, col=cols[imp_cut_index][i],border=NA)

box()


# Emulator with implausibility
plot(1, 1, xlim=c(0,1), ylim=c(-3,3), type="n", xlab="", ylab="")
upper <- em_exp + 3*sqrt(em_var); lower <- em_exp - 3*sqrt(em_var)
polygon(c(xP, rev(xP)), c(upper, rev(lower)), col=int_col, border=NA)
#lines(xP, fxP, col='red', lwd=2.5)
lines(xP, em_exp, lwd=2, col=1)#, lty=2.5)
points(xD, D, pch=21, col=1, bg=1, cex=1)

abline(h=z,lwd=2, col="green3")
abline(h=z_error,lty=2,lwd=1.5)

for(i in 1:length(xP)) rect(xP[i]-rect_wid/2,-2.9,xP[i]+rect_wid/2, -2.4, col=cols[imp_cut_index][i],border=NA)

box()
title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)
dev.off()

