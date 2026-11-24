test_that('invalid inputs are rejected and sparse inputs supported', {
  d <- fixture()
  expect_error(MSGPCA(list()),'nonempty')
  expect_error(MSGPCA(d$Y,d$z,n_shared=6,n_local=1),'dimensions')
  expect_error(MSGPCA(d$Y,d$z,lambda_shared=-1),'dimensions|lambda_shared')
  expect_error(MSGPCA(d$Y,NULL,n_shared=2,n_local=1),'Coordinates')
  a <- d$Y; colnames(a[[1]]) <- letters[1:6]
  expect_error(MSGPCA(a,d$z,n_shared=2,n_local=1),'identities')
  a <- MSGPCA(lapply(d$Y,Matrix::Matrix,sparse=TRUE),d$z,n_shared=2,n_local=1,max_iter=2,min_iter=1,verbose=FALSE)
  expect_s3_class(a,'MSGPCA')
  a <- lapply(d$Y,function(y) y*0)
  f <- MSGPCA(a,d$z,n_shared=2,n_local=4,max_iter=2,min_iter=1,verbose=FALSE)
  expect_equal(crossprod(f$Ws,f$Wl_list[[1]]),matrix(0,2,4),ignore_attr=TRUE,tolerance=1e-10)
})
test_that('seeded full fitting is reproducible and restores caller RNG', {
  set.seed(27); Y <- list(matrix(rnorm(30*120),30,120)); old <- .Random.seed
  a <- MSGPCA(Y,n_shared=2,n_local=1,lambda_shared=0,lambda_local=0,max_iter=2,min_iter=1,verbose=FALSE,seed=12)
  expect_identical(.Random.seed,old)
  b <- MSGPCA(Y,n_shared=2,n_local=1,lambda_shared=0,lambda_local=0,max_iter=2,min_iter=1,verbose=FALSE,seed=12)
  expect_identical(a$obj_history,b$obj_history); expect_identical(a$Ws,b$Ws)
})
