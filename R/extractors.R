.indices <- function(sel, ids, n, what) {
  if (is.null(sel)) return(seq_len(n))
  if (is.character(sel)) {if(is.null(ids)) stop('No ',what,' names'); out <- match(sel,ids)} else {
    if(!is.numeric(sel)||any(!is.finite(sel))||any(sel!=floor(sel))) stop('Invalid ',what,' selection')
    out <- sel
  }
  if(!length(out)||anyNA(out)||any(out<1|out>n)||anyDuplicated(out)) stop('Invalid ',what,' selection')
  out
}
.check_fit <- function(fit) if(!inherits(fit,'MSGPCA')) stop('Expected an MSGPCA fitted object')
#' Extract spatial factor scores
#' @param fit An MSGPCA fitted object.
#' @param component 'shared' or 'local'.
#' @param slices Slice indices or names; NULL selects all.
#' @param factors Factor indices or names; NULL selects all.
#' @return A slice-indexed list of score matrices. Local columns across slices
#'   need not describe the same program. No implicit score standardization.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=4,Jmix=6,Jns=4,seed=1)
#' f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=2,n_local=1,
#'             max_iter=2,min_iter=1,verbose=FALSE)
#' get_scores(f)
#' @export
get_scores <- function(fit,component='shared',slices=NULL,factors=NULL) {
  .check_fit(fit); component <- match.arg(component,c('shared','local'))
  a <- if(component=='shared') fit$Fs_list else fit$Fl_list
  si <- .indices(slices,names(a),length(a),'slice')
  lapply(a[si],function(x) x[,.indices(factors,colnames(x),ncol(x),'factor'),drop=FALSE])
}
#' Extract molecular feature loadings
#' @inheritParams get_scores
#' @return Shared loadings as one common matrix, or local loadings as a
#'   slice-indexed list. Large loadings indicate contributions, not significance.
#'   Signs and rotations are ambiguous.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=4,Jmix=6,Jns=4,seed=1)
#' f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=2,n_local=1,
#'             max_iter=2,min_iter=1,verbose=FALSE)
#' w <- get_loadings(f)
#' head(sort(abs(w[,1]),decreasing=TRUE))
#' @export
get_loadings <- function(fit,component='shared',slices=NULL,factors=NULL) {
  .check_fit(fit); component <- match.arg(component,c('shared','local'))
  if(component=='shared') return(fit$Ws[,.indices(factors,colnames(fit$Ws),ncol(fit$Ws),'factor'),drop=FALSE])
  a <- fit$Wl_list; si <- .indices(slices,names(a),length(a),'slice')
  lapply(a[si],function(x) x[,.indices(factors,colnames(x),ncol(x),'factor'),drop=FALSE])
}
