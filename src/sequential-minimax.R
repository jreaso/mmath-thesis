library(pdist)

plot_grid <- function(n){
  n = n - 1
  plot(1, type="n", xlim=c(0,1),ylim=c(0,1), xaxs="i", yaxs="i", xlab=expression("x"[1]), ylab=expression("x"[2]))
  abline(h=(0:n)/n,col="grey60")
  abline(v=(0:n)/n,col="grey60")
}

n = 25
x_grid <- seq(0, 1, len=n)
xP <- as.matrix(expand.grid("x1"=x_grid, "x2"=x_grid))

I <- function(X){
  x <- X[,1]; y <- X[,2]
  return((((x-0.4)^2 + y^2) < 0.9^2) & (y > (0.4 - 0.7*x)) & (y > (x-0.7)))
}

xP2 <- xP[I(xP),]



# Sequential minimax picking one at a time, best one each time (greedy)
sequential_minimax <- function(xP, k) {
  # Ensure k is not greater than the number of points
  if (k > nrow(xP)) {stop("k cannot be larger than the total number of points")}
  
  dist_mat <- as.matrix(dist(xP))
  
  max_dists <- apply(dist_mat, 1, max)
  xDi <- c(as.integer(which.min(max_dists))) # select first point
  
  for (j in 2:k) {
    if (j == 2) { dist_xD <- dist_mat[xDi,] } else { dist_xD <- apply(dist_mat[xDi,], 2, min) }
    nD <- nrow(dist_mat)
    dist_mat_xD <- pmin(dist_mat, matrix(dist_xD, nrow=nD, ncol=nD, byrow=TRUE))
    
    max_dists <- apply(dist_mat_xD, 1, max)
    xDi <- c(xDi, as.integer(which.min(max_dists)))
  }
  
  return(xDi)
}


xDi <- sequential_minimax(xP=xP2, k=8)
xDi
xD <- xP2[xDi,]

# plot
plot_grid(6)
points(xP2, pch=16,col="dodgerblue",cex=0.8)
points(xD, pch=16,col="red",cex=0.8)




# Choosing k best points to add but from a random sample
sequential_minimax <- function(X, k, p=3, niter=100) {
  if (k > nrow(X)) {
    stop("k cannot be larger than the total number of points")
  }
  
  XDi <- c() # Design matrix indices
  
  max_dist <- max(as.matrix(dist(X)))
  
  while(length(XDi) < k) {
    P <- min(p, k-length(XDi)) # number of points to add to XD this iteration
    
    for (i in 1:niter) {
      new_inds <- sample(setdiff(1:nrow(X), XDi), P, replace=F) # generate a sample of P random points not in XD
      XDj <- c(XDi, new_inds) # proposal for new design
      new_max_dist <- max(apply(as.matrix(pdist(X[setdiff(1:nrow(X), XDj),], X[XDj,])), 1, min))
      if(new_max_dist < max_dist){
        max_dist <- new_max_dist
        inds <- new_inds
      }
    }
    XDi <- c(XDi, inds)
  }
  return(XDi)
}


k <- 15; p <- 3
set.seed(42)
xDi <- sequential_minimax(X=xP2, k=k, p=p, niter=100000) # warning: may take a while to run
xDi
xD <- xP2[xDi,]


# plot
pdf(file="figures/2d-sequential-minimax.pdf", width=6, height=5)
par(mar=c(4,4,2,2))
plot_grid(0)
points(xP2, pch=16,col="steelblue1",cex=1)
points(xD, pch=16,col="red",cex=1.4)
dev.off()


