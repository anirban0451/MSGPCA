#' Fit Multi-Slice Graph Principal Component Analysis
#'
#' Estimates shared and slice-specific loadings and spatial scores without
#' transforming the supplied abundance matrices.
#' @param Y_list Nonempty list of observations-by-features numeric matrices;
#'   sparse matrices are supported. Common columns must have identical ordering.
#' @param locs_list List of finite coordinate matrices with matching rows, or NULL.
#' @param n_shared,n_local Nonnegative shared and local factor dimensions.
#' @param k_shared,k_local Positive nearest-neighbor counts, smaller than slice size.
#' @param lambda_shared,lambda_local Nonnegative smoothing scalars, or grids for tuning.
#' @param graph_type Either 'mutual' (default) or 'union' k-nearest-neighbor graph.
#' @param graphs_shared,graphs_local Optional lists of sparse combinatorial Laplacians.
#'   An adjacency entry is supplied as list(type='adjacency', matrix=A).
#' @param tuning Either 'none' or 'pareto'.
#' @param lambda_constraint Either 'none' or 'shared_ge_local'.
#' @param max_iter,min_iter Maximum and minimum completed update iterations.
#' @param tol Positive relative complete-objective convergence tolerance.
#' @param seed Optional nonnegative integer seed; caller RNG state is restored.
#' @param verbose Print scalar optimization progress.
#' @param covariance_cache_mb Maximum estimated MiB for persistent covariance cache.
#'   Above this budget covariances are computed one at a time; dense feature-space
#'   eigendecompositions still require quadratic feature memory.
#' @details Uses the supplied original residual-based shared-loading update and a
#'   restricted local eigendecomposition. Scores solve shifted sparse Laplacian
#'   systems. The complete summed penalized objective is monitored, including
#'   initialization (iteration zero). Coupled loading constraints mean monotonic
#'   descent is not guaranteed. Factor signs and rotations are not identifiable.
#'   Zero smoothing retains active factors; no preprocessing occurs internally.
#'   Sparse inputs may be densified for initialization, residuals and feature
#'   covariances. No observation-space inverse is constructed.
#' @return An S3 MSGPCA object with Ws, Wl_list, Fs_list, Fl_list, diagnostics,
#'   objective history, convergence and timing metadata. Tuned objects also contain
#'   optimal_lambda_shared, optimal_lambda_local and pareto_results.
#' @examples
#' sim <- simulation_gg_blocks(nside=6, Jsp=8, Jmix=12, Jns=10, seed=1)
#' fit <- MSGPCA(lapply(sim$Y_list, log1p), sim$locs_list,
#'               n_shared=2, n_local=1, max_iter=3, min_iter=1, verbose=FALSE)
#' summary(fit)
#' @export
MSGPCA <- function(Y_list, locs_list=NULL, n_shared=5L, n_local=5L,
 k_shared=4L,k_local=6L,lambda_shared=0.1,lambda_local=0.1,
 graph_type='mutual',graphs_shared=NULL,graphs_local=NULL,tuning='none',
 lambda_constraint='none',max_iter=100L,min_iter=10L,tol=1e-5,seed=NULL,
 verbose=TRUE,covariance_cache_mb=256) {
  tuning <- match.arg(tuning,c('none','pareto'))
  lambda_constraint <- match.arg(lambda_constraint,c('none','shared_ge_local'))
  args <- as.list(environment()); args$tuning <- args$lambda_constraint <- NULL
  if (tuning=='pareto') {
    args$lambda_constraint <- lambda_constraint
    ans <- do.call(tune_MSGPCA,args)
  } else {
    if (length(lambda_shared)!=1 || length(lambda_local)!=1) stop('Scalar lambdas required without tuning')
    if (lambda_constraint=='shared_ge_local' && lambda_shared<lambda_local) stop('Smoothing constraint violated')
    ans <- do.call(.fit_fixed,args)
  }
  ans$call <- match.call(); ans
}
