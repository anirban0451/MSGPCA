.valid_ids <- function(x, name) {
  if (!is.null(x) && (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))) stop('Invalid ', name)
}
.validate_data <- function(Y, locs, rs, rl, max_iter, min_iter, tol) {
  if (!is.list(Y) || !length(Y)) stop('Y_list must be a nonempty list')
  for (y in Y) {
    if (!(is.matrix(y) && is.numeric(y) || inherits(y, 'sparseMatrix') && methods::is(y, 'dMatrix')) ||
        length(dim(y)) != 2L || any(dim(y) < 1L) || any(!is.finite(y))) stop('Expression matrices must be finite numeric matrices')
    .valid_ids(colnames(y), 'feature identifiers'); .valid_ids(rownames(y), 'observation identifiers')
  }
  p <- ncol(Y[[1]])
  if (any(vapply(Y, ncol, 0L) != p)) stop('Common feature dimension required')
  cn <- lapply(Y, colnames)
  if (any(vapply(cn, function(x) !identical(x, cn[[1]]), TRUE))) stop('Feature identities and ordering must agree across slices')
  .integer(rs, 'n_shared'); .integer(rl, 'n_local')
  if (rs + rl == 0 || rs + rl > p || rs > sum(vapply(Y, nrow, 0L)) || rl > min(vapply(Y, nrow, 0L))) stop('Invalid structural factor dimensions')
  .integer(max_iter, 'max_iter', 1); .integer(min_iter, 'min_iter')
  if (min_iter > max_iter) stop('min_iter exceeds max_iter')
  .scalar(tol, 'tol'); if (tol == 0) stop('tol must be positive')
  if (!is.null(locs)) {
    if (!is.list(locs) || length(locs) != length(Y)) stop('locs_list must match slices')
    if(!is.null(names(locs)) && !is.null(names(Y)) && !identical(names(locs),names(Y))) stop('Slice ordering differs between expression and coordinates')
    dims <- vapply(locs, ncol, 0L)
    if (length(unique(dims)) != 1L || any(dims < 1)) stop('Coordinate dimensions must agree')
    for (m in seq_along(Y)) {
      z <- locs[[m]]
      if (!is.matrix(z) || !is.numeric(z) || nrow(z) != nrow(Y[[m]]) || any(!is.finite(z))) stop('Invalid coordinates')
      .valid_ids(rownames(z), 'coordinate identifiers')
      if (!is.null(rownames(z)) && !is.null(rownames(Y[[m]])) && !identical(rownames(z), rownames(Y[[m]]))) stop('Observation ordering differs')
    }
  }
  if (!is.null(names(Y))) .valid_ids(names(Y), 'slice names')
  invisible(TRUE)
}
