#' Reconstruct supplied-scale molecular abundance
#' @param fit MSGPCA fitted object.
#' @param component 'all', 'shared' or 'local'.
#' @param slices Slice names or indices, NULL for all.
#' @param features Feature names or indices, NULL for all.
#' @return Slice-indexed list of reconstructed matrices with original identifiers.
#' @details Reconstruction is on the preprocessed input scale, without inverse
#'   transformation. Requested features alone are allocated. Disabled components
#'   reconstruct to zero. No spatial extrapolation is provided.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=4,Jmix=6,Jns=4,seed=1)
#' f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=2,n_local=1,
#'             max_iter=2,min_iter=1,verbose=FALSE)
#' reconstruct_MSGPCA(f,features=1:2)
#' @export
reconstruct_MSGPCA <- function(fit,component='all',slices=NULL,features=NULL) {
  .check_fit(fit); component <- match.arg(component,c('all','shared','local'))
  si <- .indices(slices,names(fit$Fs_list),length(fit$Fs_list),'slice')
  fi <- .indices(features,rownames(fit$Ws),nrow(fit$Ws),'feature')
  ans <- lapply(si,function(m) {
    a <- tcrossprod(fit$Fs_list[[m]],fit$Ws[fi,,drop=FALSE])
    b <- tcrossprod(fit$Fl_list[[m]],fit$Wl_list[[m]][fi,,drop=FALSE])
    switch(component,all=a+b,shared=a,local=b)
  }); names(ans) <- names(fit$Fs_list)[si]; ans
}
