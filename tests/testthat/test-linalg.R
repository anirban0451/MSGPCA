test_that('permuted Cholesky solves the original A system', {
  d <- fixture(); L <- MSGPCA:::.knn_laplacian(d$z[[1]],4,'union')
  f <- MSGPCA:::.factorize(L,0.7); X <- d$Y[[1]]
  got <- MSGPCA:::.inverse_action(f,X)
  expect_equal(got,solve(diag(20)+0.7*as.matrix(L),X),tolerance=1e-10,ignore_attr=TRUE)
  expect_equal(MSGPCA:::.inverse_action(NULL,X),X)
})
test_that('null and full complement directions remain feasible', {
  W <- qr.Q(qr(matrix(1:12,6,2)))
  V <- MSGPCA:::.local_top(matrix(0,6,6),W,4)
  expect_equal(crossprod(W,V),matrix(0,2,4),tolerance=1e-12)
  expect_equal(crossprod(V),diag(4),tolerance=1e-12)
})
test_that('identical systems reuse Cholesky factors and inactive systems skip them', {
  original <- MSGPCA:::.factorize; count <- 0L
  local_mocked_bindings(.factorize=function(...) {count <<- count+1L; original(...)},.package='MSGPCA')
  d <- fixture()
  f <- MSGPCA(d$Y,d$z,n_shared=2,n_local=1,k_shared=4,k_local=4,lambda_shared=.1,lambda_local=.1,max_iter=2,min_iter=1,verbose=FALSE)
  expect_equal(count,2L)
  count <- 0L
  f <- MSGPCA(d$Y,d$z,n_shared=0,n_local=1,k_shared=4,k_local=4,max_iter=2,min_iter=1,verbose=FALSE)
  expect_equal(count,2L)
})
