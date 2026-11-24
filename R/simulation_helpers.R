.gg_geometry <- function(nside) {
  A <- matrix(0,4,36)
  A[1,c(1,6,7,8,13)+1] <- 1
  A[2,c(3,4,5,9,11,15,16,17)+1] <- 1
  A[3,c(18,24,25,30,31,32)+1] <- 1
  A[4,c(21,22,23,28,34)+1] <- 1
  F <- vapply(seq_len(4),function(l) as.vector(t(kronecker(matrix(A[l,],6,6,byrow=TRUE),matrix(1,nside/6,nside/6)))),numeric(nside^2))
  z <- as.matrix(expand.grid(x=seq_len(nside),y=seq_len(nside)))
  z[,2] <- -z[,2]
  z <- apply(z,2,function(v) (v-min(v))/(max(v)-min(v)))
  rownames(z) <- rownames(F) <- paste0('Spot',seq_len(nside^2))
  colnames(F) <- paste0('Pattern',seq_len(4))
  list(F=F,coords=z)
}
