test_that('mutual and union graphs have valid sparse Laplacians', {
  d <- fixture()
  for(type in c('mutual','union')) {
    L <- MSGPCA:::.knn_laplacian(d$z[[1]],4,type)
    expect_s4_class(L,'sparseMatrix'); expect_equal(as.matrix(L),t(as.matrix(L)))
    expect_equal(as.numeric(rowSums(L)),rep(0,20))
    expect_gte(min(eigen(as.matrix(L),symmetric=TRUE)$values),-1e-10)
    expect_true(all((L-diag(diag(L)))<=0))
  }
})
test_that('supplied adjacency and isolated vertices work', {
  A <- matrix(0,6,6); A[1,2] <- A[2,1] <- 2
  L <- MSGPCA:::.validate_graph(list(type='adjacency',matrix=A),6L)
  expect_equal(diag(L),c(2,2,0,0,0,0))
  expect_equal(MSGPCA:::.validate_graph(L,6L),L)
  expect_error(MSGPCA:::.validate_graph(diag(6),6L),'Laplacian')
  expect_error(MSGPCA:::.validate_graph(matrix(1,4,4),6L),'dimensions')
})
