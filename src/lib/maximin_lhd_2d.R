### Maximin algorithm to create 2D LHD
# More efficient to use `lhs` or `LHD` library

maximin_lhd_2d <- function(nl=16, perturb=F, niter=10000){
  # create arbitrary LHD
  x_lhd <- cbind("x1"=sample(1:nl), "x2"=sample(1:nl)) / nl - 0.5/nl
  
  # Maximin loop: performs swaps on 1st of two closest points with another random point
  for(i in 1:niter){
    dmat <- as.matrix(dist(x_lhd)) + diag(10, nl)  # creates matrix of distances between points with inflated diagonal
    closest_runs <- which(dmat==min(dmat), arr.ind=T)  # finds pairs of closest runs
    ind <- closest_runs[sample(nrow(closest_runs), 1), 1]  # chooses one of close runs at random
    swap_ind <- sample(setdiff(1:nl, ind), 1)  # randomly selects another run to swap with
    x_lhd2 <- x_lhd
    x_lhd2[ind[1], 1] <- x_lhd[swap_ind,1]
    x_lhd2[swap_ind,1] <- x_lhd[ind[1],1]  # swaps x_1 vals between 1st close run & other run
    # if min distance between points is same or better, replace LHD with the swap
    if(min(dist(x_lhd2)) >= min(dist(x_lhd))-0.00001) x_lhd <- x_lhd2
  }
  if (perturb) x_lhd <- x_lhd + matrix(runif(2*nl, min=-1/(2*nl), max=1/(2*nl)), ncol=2, nrow=nl)
  return(x_lhd)
}

