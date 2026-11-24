#' Tune MSGPCA smoothing using reconstruction and unweighted roughness
#' @param Y_list,locs_list Expression and coordinate lists as in MSGPCA.
#' @param lambda_shared,lambda_local Nonnegative finite candidate grids.
#' @param ... Fixed fitting arguments passed to the internal engine (see MSGPCA).
#' @param lambda_constraint Either 'none' or 'shared_ge_local'.
#' @param nonconverged Either 'fallback' (use finite nonconverged candidates only
#'   when none converged) or 'exclude'.
#' @return The selected MSGPCA fit, with annotated pareto_results and optimal lambdas.
#' @details Graphs are constructed once. Factorizations are scoped to one candidate,
#'   bounding memory. Candidates use identical seeded initialization; no warm starts.
#'   Converged finite fits have priority. Ineligible rows have NA normalized criteria.
#'   Disabled-component grids collapse to zero; ties follow expand.grid order.
#' @examples
#' s <- simulation_gg_blocks(nside=6, Jsp=6,Jmix=8,Jns=6,seed=1)
#' f <- tune_MSGPCA(lapply(s$Y_list,log1p),s$locs_list,c(0.1,0.3),0.1,
#'                 n_shared=2,n_local=1,max_iter=3,min_iter=1,verbose=FALSE)
#' @export
 tune_MSGPCA <- function(Y_list,locs_list=NULL,lambda_shared,lambda_local,...,
 lambda_constraint='none',nonconverged='fallback') {
  lambda_constraint <- match.arg(lambda_constraint,c('none','shared_ge_local'))
  nonconverged <- match.arg(nonconverged,c('fallback','exclude'))
  for (x in list(lambda_shared,lambda_local)) if (!is.numeric(x)||!length(x)||any(!is.finite(x)|x<0)) stop('Invalid smoothing grid')
  defaults <- formals(MSGPCA); defaults <- lapply(defaults,function(x) tryCatch(eval(x),error=function(e) NULL))
  defaults$Y_list <- Y_list; defaults$locs_list <- locs_list
  dots <- list(...)
  allowed <- setdiff(names(formals(.fit_fixed)),c('Y_list','locs_list','lambda_shared','lambda_local','prepared'))
  if (length(setdiff(names(dots),allowed))) stop('Unknown tuning argument: ',paste(setdiff(names(dots),allowed),collapse=', '))
  for (nm in names(dots)) defaults[nm] <- dots[nm]
  a <- defaults[intersect(names(defaults),names(formals(.fit_fixed)))]
  .validate_data(Y_list,locs_list,a$n_shared,a$n_local,a$max_iter,a$min_iter,a$tol)
  grid <- expand.grid(lambda_shared=if(a$n_shared) unique(lambda_shared) else 0,
                      lambda_local=if(a$n_local) unique(lambda_local) else 0)
  if (lambda_constraint=='shared_ge_local' && a$n_shared && a$n_local) grid <- grid[grid$lambda_shared>=grid$lambda_local,,drop=FALSE]
  if (!nrow(grid)) stop('No candidate satisfies restriction')
  Ls <- .graphs(Y_list,locs_list,a$graphs_shared,a$k_shared,a$graph_type,a$n_shared,max(grid$lambda_shared),TRUE)
  Ll <- if(identical(a$graphs_shared,a$graphs_local) && a$k_shared==a$k_local && a$n_shared>0 && a$n_local>0) Ls else
    .graphs(Y_list,locs_list,a$graphs_local,a$k_local,a$graph_type,a$n_local,max(grid$lambda_local),TRUE)
  prepared <- list(Ls=Ls,Ll=Ll)
  # Fixed-size selected-state retention: retain one best fit, refit winner if it changes.
  rows <- vector('list',nrow(grid))
  for(i in seq_len(nrow(grid))) {
    a$lambda_shared <- grid$lambda_shared[i]; a$lambda_local <- grid$lambda_local[i]; a$prepared <- prepared
    tm <- proc.time()[3]; error <- ''
    f <- tryCatch(do.call(.fit_fixed,a),error=function(e) {error <<- conditionMessage(e); NULL})
    h <- if(is.null(f)) rep(NA_real_,6) else unlist(tail(f$diagnostics,1)[c('reconstruction_error','shared_roughness','local_roughness','objective')])
    rows[[i]] <- data.frame(lambda_shared=a$lambda_shared,lambda_local=a$lambda_local,
      reconstruction_error=h[1],shared_roughness=h[2],local_roughness=h[3],roughness=h[2]+h[3],objective=h[4],
      converged=!is.null(f)&&f$converged,iterations=if(is.null(f)) NA_integer_ else f$iterations,
      runtime=proc.time()[3]-tm,failure=error,row.names=NULL)
  }
  d <- do.call(rbind,rows); finite <- is.finite(d$reconstruction_error)&is.finite(d$roughness)&is.finite(d$objective)&d$failure==''
  eligible <- finite & d$converged
  if (!any(eligible)&&nonconverged=='fallback') {eligible <- finite; warning('No converged candidates; selecting among finite nonconverged fits')}
  if (!any(eligible)) stop('No eligible tuning candidates')
  sel <- find_optimal_pareto_points(d[eligible,,drop=FALSE],'lambda_shared','lambda_local','reconstruction_error','roughness')
  for(nm in c('normalized_reconstruction','normalized_roughness','distance','pareto')) {
    d[[nm]] <- if(nm=='pareto') rep(NA,nrow(d)) else rep(NA_real_,nrow(d))
    d[[nm]][eligible] <- sel$results[[nm]]
  }
  d$eligible <- eligible; idx <- which(eligible)[sel$row_index]; d$selected <- seq_len(nrow(d))==idx
  a$lambda_shared <- d$lambda_shared[idx]; a$lambda_local <- d$lambda_local[idx]
  f <- do.call(.fit_fixed,a)
  f$optimal_lambda_shared <- a$lambda_shared; f$optimal_lambda_local <- a$lambda_local
  f$pareto_results <- d; f$tuning_policy <- nonconverged; f$call <- match.call(); f
}
