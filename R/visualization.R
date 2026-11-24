#' Plot MSGPCA spatial factors
#' @param fit MSGPCA fitted object.
#' @param component Shared or local scores.
#' @param slices,factors Optional slice and factor selections.
#' @param common_scale Use a common color scale across all panels.
#' @return A reusable ggplot object. Local factor indices across slices are not aligned.
#' @details Requires at least two coordinate columns; retains coordinate orientation.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=8,Jmix=12,Jns=10,seed=1)
#' f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=2,n_local=1,
#'             max_iter=2,min_iter=1,verbose=FALSE)
#' plot_factors(f)
#' @export
plot_factors <- function(fit,component='shared',slices=NULL,factors=NULL,common_scale=TRUE) {
  a <- get_scores(fit,component,slices,factors)
  if(is.null(fit$locs_list)) stop('Coordinates required')
  d <- do.call(rbind,lapply(names(a),function(m) {
    z <- fit$locs_list[[m]]; if(ncol(z)<2) stop('Two coordinate columns required')
    do.call(rbind,lapply(seq_len(ncol(a[[m]])),function(j)
      data.frame(x=z[,1],y=z[,2],value=a[[m]][,j],slice=m,factor=colnames(a[[m]])[j])))
  }))
  if(is.null(d)||!nrow(d)) stop('No factors selected')
  if(!common_scale) d$value <- ave(d$value,interaction(d$slice,d$factor),FUN=function(v) if(diff(range(v))) (v-min(v))/diff(range(v)) else rep(0,length(v)))
  ggplot2::ggplot(d,ggplot2::aes(x=x,y=y,color=value))+ggplot2::geom_point(size=0.8)+ggplot2::coord_equal()+
    ggplot2::facet_grid(slice~factor)+ggplot2::scale_color_viridis_c()+ggplot2::theme_minimal()+
    ggplot2::labs(title=paste(component,'spatial scores'),color=if(common_scale) 'Score' else 'Panel-scaled')
}
#' Plot molecular feature loadings
#' @inheritParams plot_factors
#' @param features Optional feature names or indices.
#' @param top Number of largest absolute loadings per factor when features is NULL.
#' @return Reusable ggplot object; large loadings are not significance measures.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=8,Jmix=12,Jns=10,seed=1)
#' f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=2,n_local=1,
#'             max_iter=2,min_iter=1,verbose=FALSE)
#' plot_loadings(f,top=5)
#' @export
plot_loadings <- function(fit,component='shared',slices=NULL,factors=NULL,features=NULL,top=15L) {
  .integer(top,'top',1); a <- get_loadings(fit,component,slices,factors)
  if(is.matrix(a)) a <- list(Shared=a)
  d <- do.call(rbind,lapply(names(a),function(m) {
    w <- a[[m]]
    do.call(rbind,lapply(seq_len(ncol(w)),function(j) {
      idx <- if(is.null(features)) head(order(abs(w[,j]),decreasing=TRUE),top) else .indices(features,rownames(w),nrow(w),'feature')
      data.frame(feature=if(is.null(rownames(w))) as.character(idx) else rownames(w)[idx],value=w[idx,j],factor=colnames(w)[j],slice=m)
    }))
  }))
  if(is.null(d)) stop('No factors selected')
  ggplot2::ggplot(d,ggplot2::aes(x=feature,y=value))+ggplot2::geom_col()+ggplot2::coord_flip()+
    ggplot2::facet_grid(slice~factor,scales='free_y')+ggplot2::theme_minimal()+ggplot2::labs(y='Loading')
}
#' Plot spatial cluster assignments
#' @param object An MSGPCA_clusters object.
#' @param refined Display refined labels; defaults to whether refinement was run.
#' @param slices Optional slice indices or names.
#' @return Reusable ggplot object. Local cluster labels are not aligned across slices.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=8,Jmix=12,Jns=10,seed=1)
#' f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=2,n_local=1,
#'             max_iter=2,min_iter=1,verbose=FALSE)
#' cl <- cluster_MSGPCA(f,3,method='kmeans',seed=1)
#' plot_clusters(cl)
#' @export
plot_clusters <- function(object,refined=!is.null(object$refined),slices=NULL) {
  if(!inherits(object,'MSGPCA_clusters')) stop('Expected MSGPCA_clusters')
  if(is.null(object$locs_list)) stop('Coordinates required')
  a <- if(refined) object$refined else object$raw
  if(is.null(a)) stop('Refinement was not performed')
  si <- .indices(slices,names(a),length(a),'slice')
  d <- do.call(rbind,lapply(si,function(m) {
    z <- object$locs_list[[m]]; if(ncol(z)<2) stop('Two coordinate columns required')
    lab <- a[[m]]
    if(object$parameters$component=='local') lab <- paste(names(a)[m],lab,sep=':')
    data.frame(x=z[,1],y=z[,2],cluster=lab,slice=names(a)[m])
  }))
  ggplot2::ggplot(d,ggplot2::aes(x=x,y=y,color=cluster))+ggplot2::geom_point(size=0.8)+
    ggplot2::coord_equal()+ggplot2::facet_wrap(~slice)+ggplot2::theme_minimal()
}
#' Plot Pareto candidate trade-offs
#' @param object Tuned MSGPCA object or annotated candidate data frame.
#' @return Reusable ggplot showing candidate eligibility, front and selection.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=8,Jmix=12,Jns=10,seed=1)
#' f <- tune_MSGPCA(lapply(s$Y_list,log1p),s$locs_list,c(.1,.3),.1,
#'                  n_shared=2,n_local=1,seed=1,verbose=FALSE)
#' plot_pareto(f)
#' @export
plot_pareto <- function(object) {
  d <- if(inherits(object,'MSGPCA')) object$pareto_results else object
  if(is.null(d)||!all(c('pareto','selected','reconstruction_error','roughness')%in%names(d))) stop('Annotated tuning results required')
  d$status <- ifelse(d$selected,'Selected',ifelse(is.na(d$pareto),'Ineligible',ifelse(d$pareto,'Pareto front','Candidate')))
  ggplot2::ggplot(d,ggplot2::aes(x=reconstruction_error,y=roughness,color=status))+
    ggplot2::geom_point(size=2)+ggplot2::theme_minimal()+ggplot2::labs(x='Reconstruction loss',y='Unweighted roughness',color=NULL)
}
#' Plot complete optimization objective
#' @param fit MSGPCA fitted object.
#' @return A ggplot of the complete summed objective, including iteration zero.
#' @examples
#' s <- simulation_gg_blocks(nside=6,Jsp=8,Jmix=12,Jns=10,seed=1)
#' f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=2,n_local=1,
#'             max_iter=2,min_iter=1,verbose=FALSE)
#' plot_convergence(f)
#' @export
plot_convergence <- function(fit) {
  .check_fit(fit)
  ggplot2::ggplot(fit$diagnostics,ggplot2::aes(x=iteration,y=objective))+
    ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::theme_minimal()+
    ggplot2::labs(x='Completed iteration (0 = initialization)',y='Complete penalized objective')
}
utils::globalVariables(c('x','y','value','slice','factor','feature','cluster','status','reconstruction_error','roughness','iteration','objective'))
