# Adapted from SpatialPCA, Authors L. Shang and X. Zhou, GPL >=3.
# Commit 1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66, R/SpatialPCA_utilties.R.
# Changes: validation, name-based label lookup and memory-bounded distances.
.refine_reference_neighbors <- function(location,q) {
  # One distance row at a time uses the same stats::dist arithmetic as upstream.
  # O(n^2) time, O(n*d+n*q) memory. This reference-preserving implementation
  # deliberately avoids uncertain kNN boundary ties.
  n <- nrow(location)
  t(vapply(seq_len(n),function(i) {
    d2 <- numeric(n)
    for(j in seq_len(ncol(location))) d2 <- d2+(location[i,j]-location[,j])^2
    order(sqrt(d2))[seq_len(q+1L)]
  },integer(q+1L)))
}
.refine_neighbors <- function(location,q,engine='knn') {
  if(engine=='reference') return(.refine_reference_neighbors(location,q))
  n <- nrow(location)
  nearest <- tryCatch(dbscan::kNN(location,k=q),error=function(e) NULL)
  if(is.null(nearest)||any(!is.finite(nearest$dist))) return(.refine_reference_neighbors(location,q))
  radius <- nearest$dist[,q]
  output <- matrix(0L,n,q+1L)
  # Bounded query batches: radius queries retrieve all boundary ties, not only
  # the arbitrary q neighbors returned by kNN. Exact R distances and index ties
  # select the final neighborhood. approx=0 prohibits approximate searches.
  for(start in seq.int(1L,n,by=64L)) {
    ii <- seq.int(start,min(n,start+63L))
    eps <- max(radius[ii]) + max(1,max(radius[ii]))*1e-10
    candidates <- tryCatch(dbscan::frNN(location,eps=eps,query=location[ii,,drop=FALSE],sort=FALSE,approx=0)$id,error=function(e) NULL)
    for(a in seq_along(ii)) {
      i <- ii[a]
      ids <- if(is.null(candidates)) seq_len(n) else sort(unique(c(i,candidates[[a]])))
      if(length(ids)<q+1L) ids <- seq_len(n)
      d2 <- numeric(length(ids))
      for(j in seq_len(ncol(location))) d2 <- d2+(location[i,j]-location[ids,j])^2
      output[i,] <- ids[order(sqrt(d2),ids)[seq_len(q+1L)]]
    }
  }
  output
}
#' Refine spatial labels using the SpatialPCA single-pass rule
#' @param clusterlabels Nonmissing numeric, character or factor cluster labels.
#' @param location Finite numeric observations-by-coordinate matrix, unique rows.
#' @param shape 'square' uses four neighbors; 'hexagon' uses six.
#' @param engine 'knn' uses exact radius queries to recover all boundary ties;
#'   'reference' uses rowwise distances with input-order ties.
#' @return Character labels preserving label identities and input names.
#' @details Adapted from SpatialPCA refine_cluster_10x, GPL (>=3), commit
#'   1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66. Includes the focal observation;
#'   updates sequentially once. Reassignment requires current-label frequency
#'   strictly below q/2 and maximum frequency strictly above q/2. Distance ties
#'   use input row order. Requires at least q+1 distinct observations. Shape
#'   selects neighborhood size, not physical geometry validation. Uses at most O(n^2)
#'   worst-case distance work but linear working memory rather than a dense distance
#'   matrix. The default k-NN/radius-query engine is typically faster on grids;
#'   it retrieves boundary ties, recalculates exact distances and sorts by row index.
#'   engine='reference' retains the independent rowwise implementation.
#'   Numeric/factor labels are converted to character to fix upstream positional
#'   table indexing. Missing labels and duplicate coordinates are rejected.
#' @references Shang L, Zhou X (2022). Spatially aware dimension reduction for
#'   spatial transcriptomics. Nature Communications 13, 7203.
#'   doi:10.1038/s41467-022-34879-1.
#' @examples
#' z <- as.matrix(expand.grid(x=1:3,y=1:3))
#' refine_cluster(c('a','a','a','a','b','a','a','a','a'),z)
#' @export
refine_cluster <- function(clusterlabels,location,shape='square',engine=c('knn','reference')) {
  engine <- match.arg(engine)
  shape <- match.arg(shape,c('square','hexagon')); q <- if(shape=='square') 4L else 6L
  if(!is.atomic(clusterlabels)||!is.null(dim(clusterlabels))||anyNA(clusterlabels)||!length(clusterlabels)) stop('Labels must be a nonmissing vector')
  if(!is.matrix(location)||!is.numeric(location)||nrow(location)!=length(clusterlabels)||ncol(location)<1||any(!is.finite(location))) stop('Invalid location matrix')
  if(nrow(location)<q+1L) stop('At least ',q+1L,' observations required')
  if(anyDuplicated(as.data.frame(location))) stop('Duplicate coordinates are ambiguous and unsupported')
  neighbors <- .refine_neighbors(location,q,engine)
  z <- as.character(clusterlabels); names(z) <- names(clusterlabels)
  for(i in seq_along(z)) {
    counts <- table(z[neighbors[i,]])
    if(counts[z[i]] < q/2 && max(counts)>q/2) z[i] <- names(counts)[which.max(counts)]
  }
  z
}
