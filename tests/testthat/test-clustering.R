test_that('pooled and independent clustering preserve partitions and refinement', {
  f <- fit_small()
  for(method in c('kmeans','gmm')) for(component in c('shared','local')) {
    a <- cluster_MSGPCA(f,2,method=method,component=component,seed=3)
    expect_equal(lengths(a$raw),c(A=20L,B=15L))
    for(m in 1:2) expect_identical(a$refined[[m]],refine_cluster(a$raw[[m]],f$locs_list[[m]]))
    b <- cluster_MSGPCA(f,2,method=method,component=component,refine=FALSE,seed=3)
    expect_identical(a$raw,b$raw); expect_null(b$refined)
  }
})
