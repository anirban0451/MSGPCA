test_that('shared and local only recover closed-form subspaces', {
  d <- fixture()
  for(mode in c('shared','local')) {
    rs <- if(mode=='shared') 2 else 0; rl <- if(mode=='local') 2 else 0
    f <- MSGPCA(d$Y,d$z,n_shared=rs,n_local=rl,max_iter=3,min_iter=1,verbose=FALSE)
    C <- lapply(1:2,function(m) {
      L <- MSGPCA:::.knn_laplacian(d$z[[m]],if(rs) 4 else 6,'mutual')
      crossprod(d$Y[[m]],solve(diag(nrow(L))+0.1*as.matrix(L),d$Y[[m]]))
    })
    if(rs) expect_equal(projection(f$Ws),projection(eigen(Reduce('+',C),symmetric=TRUE)$vectors[,1:2]),tolerance=1e-9,ignore_attr=TRUE)
    else for(m in 1:2) expect_equal(projection(f$Wl_list[[m]]),projection(eigen(C[[m]],symmetric=TRUE)$vectors[,1:2]),tolerance=1e-9,ignore_attr=TRUE)
  }
})
test_that('zero smoothing retains dimensions and one-component SVD', {
  d <- fixture()
  f <- MSGPCA(d$Y,NULL,n_shared=2,n_local=1,lambda_shared=0,lambda_local=0,max_iter=3,min_iter=1,verbose=FALSE)
  expect_equal(ncol(f$Fl_list[[1]]),1L); expect_equal(f$Fl_list[[1]],d$Y[[1]]%*%f$Wl_list[[1]])
  for(rs in c(0,2)) {
    f <- MSGPCA(d$Y[1],NULL,n_shared=rs,n_local=2-rs,lambda_shared=0,lambda_local=0,max_iter=3,min_iter=1,verbose=FALSE)
    V <- svd(d$Y[[1]])$v[,1:2]
    expect_equal(reconstruct_MSGPCA(f)[[1]],d$Y[[1]]%*%tcrossprod(V),tolerance=1e-9,ignore_attr=TRUE)
  }
})
