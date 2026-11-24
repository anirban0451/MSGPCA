fixture <- function() {
  set.seed(47)
  Y <- list(A=matrix(rnorm(120),20,6),B=matrix(rnorm(90),15,6))
  z <- list(matrix(runif(40),20,2),matrix(runif(30),15,2))
  list(Y=Y,z=z)
}
fit_small <- function(...) {
  d <- fixture()
  MSGPCA(d$Y,d$z,n_shared=2,n_local=1,max_iter=5,min_iter=1,verbose=FALSE,...)
}
projection <- function(W) tcrossprod(W)
upstream <- function(z,location,shape='square') {
  d <- as.matrix(dist(location)); q <- if(shape=='square') 4 else 6
  for(i in seq_along(z)) {
    a <- table(z[order(d[i,])[seq_len(q+1)]])
    if(a[z[i]]<q/2 & max(a)>q/2) z[i] <- names(a)[which.max(a)]
  }
  z
}
