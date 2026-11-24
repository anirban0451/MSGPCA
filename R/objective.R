.objective <- function(Y, Ws, Wl, Fs, Fl, Ls, Ll, ls, ll, iteration) {
  R <- Ps <- Pl <- 0
  for (m in seq_along(Y)) {
    y <- Y[[m]]; fs <- Fs[[m]]; fl <- Fl[[m]]
    nr <- sum(y*y) - 2*sum((y %*% Ws)*fs) - 2*sum((y %*% Wl[[m]])*fl) + sum(fs*fs) + sum(fl*fl)
    # Direct residual fallback when cancellation dominates.
    if (nr < 1e-9 * max(1, sum(y*y))) nr <- sum((y - tcrossprod(fs,Ws) - tcrossprod(fl,Wl[[m]]))^2)
    R <- R + nr
    Ps <- Ps + sum(fs * (Ls[[m]] %*% fs)); Pl <- Pl + sum(fl * (Ll[[m]] %*% fl))
  }
  data.frame(iteration=iteration, reconstruction_error=R, shared_roughness=Ps,
             local_roughness=Pl, weighted_shared_penalty=ls*Ps,
             weighted_local_penalty=ll*Pl, objective=R+ls*Ps+ll*Pl, relative_change=NA_real_)
}
