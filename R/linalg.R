# Internal helpers: explicit system='A' respects CHOLMOD permutations.
.with_seed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  .integer(seed, 'seed', 0)
  existed <- exists('.Random.seed', .GlobalEnv, inherits = FALSE)
  if (existed) old <- get('.Random.seed', .GlobalEnv)
  on.exit(if (existed) assign('.Random.seed', old, .GlobalEnv) else
    if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) rm('.Random.seed', envir = .GlobalEnv))
  set.seed(seed); fun()
}
.integer <- function(x, name, lower = 0) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) || x < lower)
    stop(name, ' must be an integer >= ', lower, call. = FALSE)
}
.scalar <- function(x, name, lower = 0, upper = Inf) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < lower || x > upper)
    stop(name, ' must be finite in [', lower, ', ', upper, ']', call. = FALSE)
}
.sym <- function(C) (C + t(C))/2
.project <- function(X, W) X - W %*% crossprod(W, X)
.factorize <- function(L, lambda) if (lambda == 0) NULL else
  Matrix::Cholesky(Matrix::forceSymmetric(Matrix::Diagonal(nrow(L)) + lambda * L), LDL = FALSE)
.inverse_action <- function(fac, X) {
  if (!ncol(X) || is.null(fac)) return(as.matrix(X))
  as.matrix(Matrix::solve(fac, X, system = 'A'))
}
# Complete QR provides a reliable complement even for null covariance directions.
.complement <- function(W, p) {
  if (!ncol(W)) return(diag(p))
  Q <- qr.Q(qr(W), complete = TRUE)
  Q[, seq.int(ncol(W) + 1L, p), drop = FALSE]
}
.eigen_top <- function(C, r) {
  if (!r) return(matrix(0, nrow(C), 0))
  C <- .sym(C)
  if (r < min(20L, nrow(C)/4) && nrow(C) > 100L) {
    ans <- tryCatch(irlba::irlba(C, nv = r, nu = r, tol = 1e-9), error = function(e) NULL)
    if (!is.null(ans) && max(abs(C %*% ans$u - sweep(ans$u, 2, ans$d, '*'))) < 1e-6 * max(1, max(abs(C)))) return(ans$u)
  }
  eigen(C, symmetric = TRUE)$vectors[, seq_len(r), drop = FALSE]
}
.local_top <- function(C, W, r) {
  if (!r) return(matrix(0, nrow(C), 0))
  if (!ncol(W)) return(.eigen_top(C,r))
  # Two-sided projection without constructing I-WW'.
  PC <- .project(C,W)
  PCP <- .sym(t(.project(t(PC),W)))
  V <- .project(.eigen_top(PCP,r),W)
  sv <- svd(V,nu=r,nv=0)
  if (min(sv$d)>1e-7) return(sv$u[,seq_len(r),drop=FALSE])
  # Null/repeated boundary eigenvalues can select shared-space directions.
  # Fall back to a complete orthogonal-complement basis in that case only.
  Q <- .complement(W,nrow(C))
  Q %*% .eigen_top(crossprod(Q,C %*% Q),r)
}
.svd_right <- function(Y, r) {
  if (!r) return(matrix(0, ncol(Y), 0))
  if (r < min(dim(Y))/4 && min(dim(Y)) > 40) {
    a <- tryCatch(irlba::irlba(Y, nv = r, nu = 0, tol = 1e-9), error = function(e) NULL)
    if (!is.null(a)) return(a$v)
  }
  svd(Y, nu = 0, nv = r)$v[, seq_len(r), drop = FALSE]
}
