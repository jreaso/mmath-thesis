source("src/lib/simple_BL_emulator.R")

# Non-stationary function
f <- function(x) {
  x + sin(3*pi*x)/3
}
  
xP <- seq(-3.25, 3.25, len=1001)
fxP <- lapply(xP, f)

xD <- seq(-1.5, 1.5, l=10)
D <- f(xD)


# Constant Prior Expectation
em_out <- simple_BL_emulator(xD=as.matrix(xD, ncol=1), D=D, xP=as.matrix(xP, ncol=1),
                             theta=.35, sig=.5, nugget=0, mu=0)
em_exp <- as.vector(em_out[["ExpD_f(x)"]])
em_var <- as.vector(em_out[["VarD_f(x)"]])


int_col <- rgb(.9,.9,.9)

pdf(file="figures/extrap-probs.pdf", width=6, height=5)
  plot(1, 1, xlim=c(-3,3), ylim=c(-3,3), type="n", xlab="", ylab="")
  upper <- em_exp + 3*sqrt(em_var); lower <- em_exp - 3*sqrt(em_var)
  polygon(c(xP, rev(xP)), c(upper, rev(lower)), col=int_col, border=NA)
  lines(xP, fxP, col='red', lwd=2.5)
  lines(xP, em_exp, col=1, lwd=2, lty=2.5)
  points(xD, D, pch=21, col=1, bg=1, cex=1)
  box()
  title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)


em_out <- simple_BL_emulator(xD=as.matrix(xD, ncol=1), D=D, xP=as.matrix(xP, ncol=1),
                             theta=.35, sig=.5, nugget=0, mu=function(x) x)
em_exp <- as.vector(em_out[["ExpD_f(x)"]])
em_var <- as.vector(em_out[["VarD_f(x)"]])


int_col <- rgb(.9,.9,.9)

plot(1, 1, xlim=c(-3,3), ylim=c(-3,3), type="n", xlab="", ylab="")
upper <- em_exp + 3*sqrt(em_var); lower <- em_exp - 3*sqrt(em_var)
polygon(c(xP, rev(xP)), c(upper, rev(lower)), col=int_col, border=NA)
lines(xP, fxP, col='red', lwd=2.5)
lines(xP, em_exp, col=1, lwd=2, lty=2.5)
points(xD, D, pch=21, col=1, bg=1, cex=1)
box()
title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)
dev.off()


