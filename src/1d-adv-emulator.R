source("src/lib/simple_BL_emulator.R")
source("src/lib/advanced_BL_emulator.R")

# Non-stationary function
f <- function(x) {
  x^2 / 2 + sin(3*pi*x)/5 * (1+x)
}

xP <- seq(-3.25, 3.25, len=1001)
fxP <- lapply(xP, f)

xD <- c(-3, -2.7, -2.4, -2.1, -1.8, -1.5, 0.3, 0.6, 0.9, 1.2, 1.5, 2.5, 2.8)
D <- f(xD)

# Simple Emulator
em_out <- simple_BL_emulator(xD=as.matrix(xD, ncol=1), D=D, xP=as.matrix(xP, ncol=1),
                             theta=.35, sig=sqrt(2), nugget=0, mu=0)
em_exp <- as.vector(em_out[["ExpD_f(x)"]])
em_var <- as.vector(em_out[["VarD_f(x)"]])

int_col <- rgb(.9,.9,.9)

pdf(file="figures/1d-adv-emulator.pdf", width=6, height=5)
plot(1, 1, xlim=c(-3,3), ylim=c(-1,5), type="n", xlab="", ylab="")
upper <- em_exp + 3*sqrt(em_var); lower <- em_exp - 3*sqrt(em_var)
polygon(c(xP, rev(xP)), c(upper, rev(lower)), col=int_col, border=NA)
lines(xP, fxP, col='red', lwd=2.5)
lines(xP, em_exp, col=1, lwd=2, lty=2.5)
points(xD, D, pch=21, col=1, bg=1, cex=1)
box()
title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)


### Applying advanced emulator

g <- list(function(x) 1, function(x) x, function(x) x^2)
mu_beta <- rep(0, length(g))
Sigma_beta <- matrix(rep(1, length(g)^2), nrow=length(g)) + 7*diag(length(g))

em_out <- advanced_BL_emulator(xD=as.matrix(xD, ncol=1), D=D, xP=as.matrix(xP, ncol=1),
                               mu_beta=mu_beta , g=g, Sigma_beta=Sigma_beta,
                               theta=.3, sig=0.4, nugget=0)

em_exp <- as.vector(em_out[["ExpD_f(x)"]])
em_var <- as.vector(diag(em_out[["VarD_f(x)"]]))
ED_beta <- as.vector(em_out[["ExpD_beta"]])
GP <- as.matrix(em_out[["GP"]])



plot(1, 1, xlim=c(-3,3), ylim=c(-1,5), type="n", xlab="", ylab="")
upper <- em_exp + 3*sqrt(em_var); lower <- em_exp - 3*sqrt(em_var)
polygon(c(xP, rev(xP)), c(upper, rev(lower)), col=int_col, border=NA)
lines(xP, as.vector(GP %*% ED_beta), col="royalblue2", lwd=2, lty=2.5)
lines(xP, fxP, col='red', lwd=2.5)
lines(xP, em_exp, col=1, lwd=2, lty=2.5)
points(xD, D, pch=21, col=1, bg=1, cex=1)
box()
title(ylab="f(x)", xlab="x", line=2.5, cex.lab=1)
dev.off()

