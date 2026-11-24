#' Detect spatial domains from MSGPCA scores
#' @param fit MSGPCA fitted object.
#' @param n_clusters Positive requested number of clusters.
#' @param method 'gmm' (mclust covariance-model selection) or 'kmeans'.
#' @param component 'shared' pooled across slices, or 'local' clustered separately.
#' @param factors Factor indices to use, NULL for all. No truth-based selection.
#' @param refine Apply the single-pass SpatialPCA-derived refinement.
#' @param refine_shape 'square' or 'hexagon', one value or one per slice.
#' @param nstart K-means random starts (default 25).
#' @param seed Optional seed with caller RNG restored.
#' @return MSGPCA_clusters list with raw and refined slice label lists (refined
#'   is NULL when disabled), selected GMM models, parameters, coordinates,
#'   reassignment counts and timings. Local cluster identities are not aligned.
#' @details Shared scores are clustered jointly; local coordinates are never
#'   pooled. No implicit standardization. GMM failures raise informative errors;
#'   there is no silent method fallback. Both algorithms use refine_cluster.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=8,Jmix=12,Jns=10,seed=1)
#' f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=2,n_local=1,
#'             max_iter=2,min_iter=1,verbose=FALSE)
#' cl <- cluster_MSGPCA(f,3,method='kmeans',seed=1)
#' @export
cluster_MSGPCA <- function(fit,n_clusters,method=c('gmm','kmeans'),component='shared',
 factors=NULL,refine=TRUE,refine_shape='square',nstart=25L,seed=NULL) {
  .check_fit(fit); method <- match.arg(method); component <- match.arg(component,c('shared','local'))
  .integer(n_clusters,'n_clusters',1); .integer(nstart,'nstart',1)
  if(!is.logical(refine)||length(refine)!=1||is.na(refine)) stop('refine must be logical')
  X <- get_scores(fit,component,factors=factors); M <- length(X)
  if(any(vapply(X,ncol,0L)==0)) stop('No active selected factors')
  if(refine && is.null(fit$locs_list)) stop('Coordinates required for refinement')
  if(!length(refine_shape)%in%c(1L,M)||any(!refine_shape%in%c('square','hexagon'))) stop('Invalid refinement shapes')
  shapes <- rep(refine_shape,length.out=M)
  cluster <- function(x) {
    if(n_clusters>nrow(x)) stop('More clusters than observations')
    if(method=='kmeans') return(list(labels=stats::kmeans(x,centers=n_clusters,nstart=nstart)$cluster,model=NA_character_))
    f <- mclust::Mclust(x,G=n_clusters,verbose=FALSE)
    if(is.null(f)||is.null(f$classification)||anyNA(f$classification)||length(unique(f$classification))!=n_clusters) stop('GMM fitting failed to produce requested clusters')
    list(labels=f$classification,model=f$modelName)
  }
  tm <- proc.time()[3]
  result <- .with_seed(seed,function() {
    if(component=='shared') {
      a <- cluster(do.call(rbind,X)); end <- cumsum(vapply(X,nrow,0L)); begin <- c(0L,head(end,-1))
      list(raw=lapply(seq_len(M),function(m) a$labels[seq.int(begin[m]+1,end[m])]),models=a$model)
    } else {
      a <- lapply(X,cluster); list(raw=lapply(a,`[[`,'labels'),models=vapply(a,`[[`,'','model'))
    }
  })
  raw <- lapply(result$raw,as.character); names(raw) <- names(X)
  ct <- proc.time()[3]-tm; tm <- proc.time()[3]
  refined <- if(refine) lapply(seq_len(M),function(m) refine_cluster(raw[[m]],fit$locs_list[[m]],shapes[m])) else NULL
  if(refine) names(refined) <- names(X)
  changed <- if(refine) vapply(seq_len(M),function(m) sum(raw[[m]]!=refined[[m]]),0L) else rep(0L,M)
  structure(list(raw=raw,refined=refined,locs_list=fit$locs_list,slice_ids=names(X),
    factors=colnames(X[[1]]),parameters=list(method=method,component=component,n_clusters=n_clusters,refine=refine,nstart=nstart,seed=seed),
    gmm_model=result$models,refine_shape=shapes,reassigned=changed,
    reassigned_proportion=changed/vapply(X,nrow,0L),clustering_runtime=ct,refinement_runtime=proc.time()[3]-tm),class='MSGPCA_clusters')
}
