test_that('components sum and selections preserve dimensions and names', {
  f <- fit_small(); a <- reconstruct_MSGPCA(f); s <- reconstruct_MSGPCA(f,'shared'); l <- reconstruct_MSGPCA(f,'local')
  for(m in 1:2) expect_equal(a[[m]],s[[m]]+l[[m]])
  expect_equal(reconstruct_MSGPCA(f,slices='B',features=c(2,4))[[1]],a[[2]][,c(2,4),drop=FALSE])
  expect_equal(fitted(f),a)
  expect_s3_class(plot_convergence(f),'ggplot')
  expect_s3_class(plot_factors(f),'ggplot'); expect_s3_class(plot_loadings(f),'ggplot')
})
