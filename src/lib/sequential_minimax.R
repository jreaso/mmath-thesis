sequential_minimax <- function(X, k, p=3, niter=100, XDi=c()) {
  if (k > nrow(X)) {
    stop("k cannot be larger than the total number of points")
  }
  
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
