.fit_fixed <- function(Y_list, locs_list, n_shared, n_local, k_shared, k_local,
 lambda_shared, lambda_local, graph_type, graphs_shared, graphs_local,
 max_iter, min_iter, tol, seed, verbose, covariance_cache_mb=256, prepared=NULL) {
  .with_seed(seed,function() {
  .validate_data(Y_list,locs_list,n_shared,n_local,max_iter,min_iter,tol)
  .scalar(lambda_shared,'lambda_shared'); .scalar(lambda_local,'lambda_local')
  .scalar(covariance_cache_mb,'covariance_cache_mb')
  .integer(k_shared,'k_shared',1); .integer(k_local,'k_local',1)
  graph_type <- match.arg(graph_type,c('mutual','union'))
  Y <- Y_list; M <- length(Y); p <- ncol(Y[[1]])
  sn <- names(Y); if (is.null(sn)) sn <- paste0('Slice',seq_len(M))
  names(Y) <- sn; if (!is.null(locs_list)) names(locs_list) <- sn
  start <- proc.time()[3]
  if (is.null(prepared)) {
    Ls <- .graphs(Y,locs_list,graphs_shared,k_shared,graph_type,n_shared,lambda_shared)
    Ll <- if (identical(graphs_shared,graphs_local) && k_shared==k_local && n_shared>0 && n_local>0) Ls else
      .graphs(Y,locs_list,graphs_local,k_local,graph_type,n_local,lambda_local)
  } else { Ls <- prepared$Ls; Ll <- prepared$Ll }
  graph_time <- proc.time()[3]-start; t0 <- proc.time()[3]
  Ks <- if(n_shared) lapply(Ls,.factorize,lambda=lambda_shared) else rep(list(NULL),M)
  Kl <- if (identical(Ls,Ll) && lambda_shared==lambda_local) Ks else if(n_local) lapply(Ll,.factorize,lambda=lambda_local) else rep(list(NULL),M)
  factor_time <- proc.time()[3]-t0
  # At most two p-by-p covariances per slice. Stream the shared base sum.
  cache <- (M+1)*p*p*8 <= covariance_cache_mb*1024^2
  Cs <- if(cache && n_shared) matrix(0,p,p) else NULL; Cl <- vector('list',M)
  t0 <- proc.time()[3]
  if (cache) for (m in seq_len(M)) {
    if (n_shared) Cs <- Cs + as.matrix(crossprod(Y[[m]],.inverse_action(Ks[[m]],Y[[m]])))
    if (n_local) Cl[[m]] <- as.matrix(crossprod(Y[[m]],.inverse_action(Kl[[m]],Y[[m]])))
  }
  precompute_time <- proc.time()[3]-t0; t0 <- proc.time()[3]
  state <- .with_seed(seed, function() {
    Ws <- .svd_right(do.call(rbind,Y),n_shared)
    Wl <- lapply(Y,function(y) {
      if (!n_local) return(matrix(0,p,0))
      yr <- y - tcrossprod(y %*% Ws,Ws)
      # Residual covariance restricted to a reliable complement (SVD-equivalent).
      .local_top(as.matrix(crossprod(yr)),Ws,n_local)
    })
    list(Ws=Ws,Wl=Wl)
  })
  Ws <- state$Ws; Wl <- state$Wl
  score <- function(Ws,Wl) list(Fs=lapply(seq_len(M),function(m) .inverse_action(Ks[[m]],Y[[m]] %*% Ws)),
                              Fl=lapply(seq_len(M),function(m) .inverse_action(Kl[[m]],Y[[m]] %*% Wl[[m]])))
  s <- score(Ws,Wl); Fs <- s$Fs; Fl <- s$Fl
  history <- .objective(Y,Ws,Wl,Fs,Fl,Ls,Ll,lambda_shared,lambda_local,0L)
  if(!is.finite(history$objective)) stop('Nonfinite initialization objective')
  init_time <- proc.time()[3]-t0; t0 <- proc.time()[3]; converged <- FALSE
  for (iter in seq_len(max_iter)) {
    if (n_shared) {
      C <- if (cache) Cs else matrix(0,p,p)
      for (m in seq_len(M)) {
        if (cache) {
          if (n_local) {
            KF <- .inverse_action(Ks[[m]],Fl[[m]])
            B <- as.matrix(crossprod(Y[[m]],KF)); D <- crossprod(Fl[[m]],KF)
            C <- C - tcrossprod(B,Wl[[m]]) - tcrossprod(Wl[[m]],B) + Wl[[m]] %*% D %*% t(Wl[[m]])
          }
        } else {
          R <- Y[[m]] - tcrossprod(Fl[[m]],Wl[[m]])
          C <- C + as.matrix(crossprod(R,.inverse_action(Ks[[m]],R)))
        }
      }
      Ws <- .eigen_top(C,n_shared)
    }
    if (n_local) for (m in seq_len(M)) {
      C <- if (cache) Cl[[m]] else as.matrix(crossprod(Y[[m]],.inverse_action(Kl[[m]],Y[[m]])))
      Wl[[m]] <- .local_top(C,Ws,n_local)
    }
    s <- score(Ws,Wl); Fs <- s$Fs; Fl <- s$Fl
    h <- .objective(Y,Ws,Wl,Fs,Fl,Ls,Ll,lambda_shared,lambda_local,iter)
    if(!is.finite(h$objective)) stop('Nonfinite optimization objective')
    h$relative_change <- abs(h$objective-tail(history$objective,1))/max(abs(tail(history$objective,1)),.Machine$double.eps)
    history <- rbind(history,h)
    if (verbose) message(sprintf('Iteration %d | objective %.8g | relative change %.3g',iter,h$objective,h$relative_change))
    if (iter >= min_iter && h$relative_change < tol) {converged <- TRUE; break}
  }
  for (m in seq_len(M)) {
    rownames(Fs[[m]]) <- rownames(Fl[[m]]) <- rownames(Y[[m]])
    rownames(Wl[[m]]) <- colnames(Y[[m]])
    colnames(Fs[[m]]) <- colnames(Ws) <- if (n_shared) paste0('Shared',seq_len(n_shared)) else character()
    colnames(Fl[[m]]) <- colnames(Wl[[m]]) <- if (n_local) paste0('Local',seq_len(n_local)) else character()
  }
  rownames(Ws) <- colnames(Y[[1]])
  names(Fs) <- names(Fl) <- names(Wl) <- sn
  increases <- diff(history$objective); material <- increases > 1e-8*pmax(1,abs(head(history$objective,-1)))
  structure(list(Ws=Ws,Wl_list=Wl,Fs_list=Fs,Fl_list=Fl,locs_list=locs_list,
    diagnostics=history,obj_history=history$objective,converged=converged,iterations=iter,
    relative_change=tail(history$relative_change,1),termination_reason=if(converged) 'relative_objective_change' else 'maximum_iterations',
    objective_increases=list(number=sum(material),maximum=if(any(material)) max(increases[material]) else 0, magnitudes=increases[material]),
    graph_diagnostics=list(shared=lapply(Ls,function(L) list(isolated=sum(Matrix::diag(L)==0),edges=sum(L<0)/2)),local=lapply(Ll,function(L) list(isolated=sum(Matrix::diag(L)==0),edges=sum(L<0)/2))),
    timings=c(graph=graph_time,factorization=factor_time,precomputation=precompute_time,initialization=init_time,optimization=proc.time()[3]-t0,total=proc.time()[3]-start),
    parameters=list(n_shared=n_shared,n_local=n_local,k_shared=k_shared,k_local=k_local,lambda_shared=lambda_shared,lambda_local=lambda_local,graph_type=graph_type,max_iter=max_iter,min_iter=min_iter,tol=tol,seed=seed,covariance_cache=cache,covariance_estimate_bytes=(M+1)*p*p*8),call=NULL),class='MSGPCA')
  })
}
