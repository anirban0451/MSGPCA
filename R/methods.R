#' @export
#' @method print MSGPCA
print.MSGPCA <- function(x,...) {
  cat('MSGPCA 0.1.0:',length(x$Fs_list),'slices;',ncol(x$Ws),'shared and',ncol(x$Wl_list[[1]]),'local factors\n')
  cat(x$termination_reason,'after',x$iterations,'iterations; objective',tail(x$obj_history,1),'\n')
  invisible(x)
}
#' @export
#' @method summary MSGPCA
summary.MSGPCA <- function(object,...) {
  list(parameters=object$parameters,converged=object$converged,termination_reason=object$termination_reason,
       iterations=object$iterations,final_diagnostics=tail(object$diagnostics,1),
       objective_increases=object$objective_increases,timings=object$timings)
}
#' @export
fitted.MSGPCA <- function(object,...) reconstruct_MSGPCA(object,...)
#' @export
plot.MSGPCA <- function(x,y=NULL,...) plot_convergence(x,...)
