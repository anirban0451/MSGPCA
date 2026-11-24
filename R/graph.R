.validate_graph <- function(g, n) {
  kind <- if (is.list(g)) g$type else 'laplacian'
  if (is.list(g)) g <- g$matrix
  if (!kind %in% c('adjacency', 'laplacian')) stop('Graph type must be adjacency or laplacian')
  if (!(is.matrix(g) || inherits(g, 'Matrix')) || !identical(dim(g), c(n,n)) || any(!is.finite(g))) stop('Invalid graph dimensions or values')
  g <- methods::as(Matrix::Matrix(g, sparse = TRUE), 'generalMatrix')
  if (max(abs(g - t(g))) > 1e-9) stop('Graph must be symmetric')
  if (kind == 'adjacency') {
    if (any(g < 0) || any(Matrix::diag(g) != 0)) stop('Adjacency weights must be nonnegative with zero diagonal')
    g <- Matrix::Diagonal(x = Matrix::rowSums(g)) - g
  } else {
    off <- g - Matrix::Diagonal(x = Matrix::diag(g))
    if (any(Matrix::diag(g) < -1e-9) || any(off > 1e-9) || max(abs(Matrix::rowSums(g))) > 1e-8 * max(1,max(abs(g)))) stop('Invalid combinatorial Laplacian')
  }
  Matrix::forceSymmetric(g)
}
.knn_laplacian <- function(x, k, type) {
  .integer(k, 'neighborhood size', 1)
  n <- nrow(x); if (k >= n) stop('Neighborhood size must be less than observations')
  nn <- dbscan::kNN(x, k = k)$id
  A <- Matrix::sparseMatrix(i = rep(seq_len(n), each=k), j = as.vector(t(nn)), x=1, dims=c(n,n))
  A <- if (type == 'mutual') A * t(A) else (A + t(A)) > 0
  A <- A * 1
  Matrix::forceSymmetric(Matrix::Diagonal(x=Matrix::rowSums(A)) - A)
}
.graphs <- function(Y, locs, supplied, k, type, active, lambda, force = FALSE) {
  M <- length(Y)
  if (!is.null(supplied) && (!is.list(supplied) || length(supplied) != M)) stop('Graphs must be a list matching slices')
  lapply(seq_len(M), function(m) {
    if (!active) return(Matrix::Diagonal(nrow(Y[[m]]), x=0))
    if (!is.null(supplied)) return(.validate_graph(supplied[[m]], nrow(Y[[m]])))
    if (is.null(locs)) {
      if (lambda > 0 || force) stop('Coordinates or supplied graphs required')
      return(Matrix::Diagonal(nrow(Y[[m]]), x=0))
    }
    .knn_laplacian(locs[[m]], k, type)
  })
}
