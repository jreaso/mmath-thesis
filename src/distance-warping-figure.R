library(plot3D)

v <- function(X) {
  if (is.matrix(X)) {
    x <- X[,1]; y <- X[,2]
  } else {
    x <- X[1]; y <- X[2]
  }
  (x > 0)*x^2
}

x_mesh <- seq(-0.1,1, 0.1)
y_mesh <- seq(-0.5,0.2, 0.1)
mesh <- as.matrix(expand.grid("x"=x_mesh, "y"=y_mesh))

z <- matrix(v(mesh), nrow=length(x_mesh), ncol=length(y_mesh))


# Simple diagram to show why the distances warp
pdf(file="figures/dist-warp.pdf", width=6, height=4.5)
  par(mar=c(0,0,0,0))
  persp3D(x_mesh, y_mesh, z,
          theta=-140, phi=10, col="grey90", axes=FALSE, box=FALSE, lighting=TRUE, facets=FALSE, expand=0.5, alpha=0.5)

  lines3D(c(0.2, 0.2), c(0,0), c(0.2^2, 0.8^2), col = "grey30", add=TRUE, lwd=3, lty=2)
  lines3D(c(0.2, 0.8), c(0,0), c(0.8^2, 0.8^2), col = "dodgerblue3", add=TRUE, lwd=4)
  lines(x=c(-0.2, 0.023), y=c(0.04, -0.18), col="red3", lwd=4)
  points(x=-0.2, y=0.04, col="grey15", bg="grey15", cex=3, pch=21)
  points(x=0.023, y=-0.18, col="grey15", bg="grey15", cex=3, pch=21)
  points(x=0.023, y=0.01, cex=3, col="grey15", bg="grey15", pch=21)
dev.off()
