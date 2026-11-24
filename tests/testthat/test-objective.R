test_that('both covariance identities agree with direct residual algebra', {
  d <- fixture(); f <- fit_small(); Y <- d$Y[[1]]; Ws <- f$Ws; Wl <- f$Wl_list[[1]]
  F <- f$Fl_list[[1]]; Fs <- f$Fs_list[[1]]
  L <- MSGPCA:::.knn_laplacian(d$z[[1]],4,'mutual'); K <- solve(diag(nrow(Y))+0.2*as.matrix(L))
  P <- diag(ncol(Y))-tcrossprod(Ws); C <- crossprod(Y,K%*%Y)
  R <- Y-tcrossprod(Fs,Ws)
  expect_equal(P%*%crossprod(R,K%*%R)%*%P,P%*%C%*%P,tolerance=1e-10,ignore_attr=TRUE)
  B <- crossprod(Y,K%*%F); D <- crossprod(F,K%*%F); R <- Y-tcrossprod(F,Wl)
  expect_equal(C-tcrossprod(B,Wl)-tcrossprod(Wl,B)+Wl%*%D%*%t(Wl),crossprod(R,K%*%R),tolerance=1e-10,ignore_attr=TRUE)
  direct <- sum((Y-tcrossprod(Fs,Ws)-tcrossprod(F,Wl))^2)
  fast <- sum(Y^2)-2*sum((Y%*%Ws)*Fs)-2*sum((Y%*%Wl)*F)+sum(Fs^2)+sum(F^2)
  expect_equal(direct,fast,tolerance=1e-11)
})
