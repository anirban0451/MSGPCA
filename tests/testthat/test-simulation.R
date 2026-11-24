test_that('geometry, RNG, labels and metadata are consistent', {
  set.seed(42); old <- .Random.seed
  a <- simulation_gg_blocks(nside=12,seed=101)
  expect_identical(old,.Random.seed)
  expect_identical(a,simulation_gg_blocks(nside=12,seed=101))
  expect_equal(dim(a$Y_list[[1]]),c(144L,400L))
  expect_true(all(a$Y_list[[1]]>=0 & a$Y_list[[1]]==floor(a$Y_list[[1]])))
  expect_equal(lengths(a$true_cluster_labels),c(Slice1=144L,Slice2=144L))
  expect_true(all(rowSums(a$true_shared_factors)<=1))
  expect_equal(length(a$local_gene_indices[[1]]),10L)
  z <- a$locs_list[[1]]
  expect_identical(a$perturbations$masks$diagonal,z[,1]+z[,2]<=0.9)
  expect_true(all(a$true_local_patterns[[1]]==0))
  expect_equal(a$true_local_patterns[[2]][,1],as.numeric(a$perturbations$masks[[1]]),ignore_attr=TRUE)
})
test_that('scenarios preserve unaffected slices and category boundaries', {
  base <- simulation_gg_blocks(nside=6,n_slices=5,scenario='none',seed=2)
  for(scenario in c('slab','wave','centroid_shift')) {
    s <- simulation_gg_blocks(nside=6,n_slices=5,scenario=scenario,seed=2)
    expect_identical(s$Y_list[[1]],base$Y_list[[1]])
    if(scenario=='centroid_shift') {
      expect_identical(s$Y_list[[2]],base$Y_list[[2]])
      expect_true(all(s$Y_list[[3]]>=base$Y_list[[3]]))
      expect_identical(s$Y_list[[3]][,-(1:50)],base$Y_list[[3]][,-(1:50)])
    } else {
      untouched <- setdiff(1:400,unlist(s$local_gene_indices))
      expect_identical(s$Y_list[[2]][,untouched],base$Y_list[[2]][,untouched])
    }
  }
  expect_error(simulation_gg_blocks(n_slices=2,scenario='centroid_shift'),'five|5')
  expect_error(simulation_gg_blocks(nside=7),'six')
  expect_error(simulation_gg_blocks(Lns=0),'composition')
  s <- simulation_gg_blocks(nside=6,Jsp=3,Jmix=0,Jns=0,Lns=0,scenario='none')
  expect_equal(dim(s$Y_list[[1]]),c(36L,3L))
  s <- simulation_gg_blocks(nside=6,Jsp=0,Jmix=0,Jns=5,mix_frac_red=0,mix_frac_green=0)
  expect_length(s$local_gene_indices[[1]],0L)
})
