test_that('normalization, dominance and ties match hand calculations', {
  d <- data.frame(s=1:4,l=1:4,R=c(1,2,4,4),P=c(4,2,1,4))
  a <- find_optimal_pareto_points(d,'s','l','R','P')
  expect_equal(a$row_index,2L)
  expect_equal(a$results$normalized_reconstruction,c(0,1/3,1,1))
  expect_equal(a$results$pareto,c(TRUE,TRUE,TRUE,FALSE))
  d$R <- d$P <- 1
  a <- find_optimal_pareto_points(d,'s','l','R','P')
  expect_equal(a$row_index,1L); expect_true(all(a$results$distance==0))
  expect_error(find_optimal_pareto_points(d[0,],'s','l','R','P'),'Invalid')
})
test_that('tuning and selected model agree including disabled grids', {
  d <- fixture()
  f <- suppressWarnings(tune_MSGPCA(d$Y,d$z,c(0.1,0.3),c(0.1,0.2),n_shared=2,n_local=1,max_iter=3,min_iter=1,verbose=FALSE,seed=1,lambda_constraint='shared_ge_local'))
  expect_true(all(f$pareto_results$lambda_shared>=f$pareto_results$lambda_local))
  s <- f$pareto_results[f$pareto_results$selected,]
  expect_equal(f$parameters$lambda_shared,s$lambda_shared)
  expect_equal(tail(f$obj_history,1),s$objective,tolerance=1e-9)
  f <- tune_MSGPCA(d$Y,d$z,c(0.1,0.3),c(0.1,0.2),n_shared=2,n_local=0,max_iter=3,min_iter=1,verbose=FALSE)
  expect_equal(nrow(f$pareto_results),2L)
  expect_error(tune_MSGPCA(d$Y,d$z,0.1,0.2,n_shared=2,n_local=1,max_iter=1,min_iter=1,tol=1e-15,verbose=FALSE,nonconverged='exclude'),'No eligible')
})
