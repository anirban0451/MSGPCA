test_that('final state, constraints, equations and objective agree independently', {
  d <- fixture(); f <- fit_small(lambda_shared=0.3,lambda_local=0.1)
  expect_equal(unname(crossprod(f$Ws)),diag(2),tolerance=1e-10)
  O <- 0
  for(m in 1:2) {
    Ls <- MSGPCA:::.knn_laplacian(d$z[[m]],4,'mutual')
    Ll <- MSGPCA:::.knn_laplacian(d$z[[m]],6,'mutual')
    expect_equal(unname(crossprod(f$Wl_list[[m]])),diag(1),tolerance=1e-10)
    expect_equal(unname(crossprod(f$Ws,f$Wl_list[[m]])),matrix(0,2,1),tolerance=1e-10)
    expect_equal(as.matrix((diag(nrow(Ls))+0.3*Ls)%*%f$Fs_list[[m]]),d$Y[[m]]%*%f$Ws,tolerance=1e-9)
    expect_equal(as.matrix((diag(nrow(Ll))+0.1*Ll)%*%f$Fl_list[[m]]),d$Y[[m]]%*%f$Wl_list[[m]],tolerance=1e-9)
    O <- O+sum((d$Y[[m]]-reconstruct_MSGPCA(f)[[m]])^2)+0.3*sum(f$Fs_list[[m]]*(Ls%*%f$Fs_list[[m]]))+0.1*sum(f$Fl_list[[m]]*(Ll%*%f$Fl_list[[m]]))
  }
  expect_equal(tail(f$obj_history,1),O,tolerance=1e-10)
  inc <- diff(f$obj_history); mat <- inc>1e-8*pmax(1,abs(head(f$obj_history,-1)))
  expect_equal(f$objective_increases$number,sum(mat)); expect_equal(f$objective_increases$magnitudes,inc[mat])
  expect_equal(nrow(f$diagnostics),f$iterations+1L)
})
test_that('maximum and convergence termination states are coherent', {
  d <- fixture()
  f <- MSGPCA(d$Y,d$z,n_shared=2,n_local=1,max_iter=1,min_iter=1,tol=1e-15,verbose=FALSE)
  expect_equal(f$iterations,1L); expect_false(f$converged)
  f <- MSGPCA(d$Y,d$z,n_shared=2,n_local=0,max_iter=4,min_iter=2,verbose=FALSE)
  expect_true(f$converged); expect_equal(f$iterations,2L)
})
test_that('covariance-cache and streamed residual engines agree', {
  a <- fit_small(covariance_cache_mb=256,seed=1)
  b <- fit_small(covariance_cache_mb=0,seed=1)
  expect_equal(projection(a$Ws),projection(b$Ws),tolerance=1e-8,ignore_attr=TRUE)
  expect_equal(a$obj_history,b$obj_history,tolerance=1e-8)
})
