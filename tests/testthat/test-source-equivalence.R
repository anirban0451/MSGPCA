test_that('simulation baseline and wave match unchanged upstream RNG sequence', {
  # Functions copied verbatim from supplied source into a test fixture.
  env <- new.env(parent=baseenv())
  # Stats helpers referenced by original unqualified functions.
  for(nm in c('rbinom','rnbinom','model.matrix')) assign(nm,get(nm,asNamespace('stats')),env)
  source(test_path('fixtures','original_simulation.R'),local=env)
  set.seed(101); baseline <- env$simulate_dataset()
  orig <- env$inject_targeted_local_wave_features(baseline$Ylist,baseline$locslist,baseline$true_spfac,intensity_mean=8)
  s <- simulation_gg_blocks(seed=101)
  expect_identical(unname(s$Y_list),orig$Ylist)
  expect_identical(unname(s$locs_list),baseline$locslist)
  expect_equal(unname(s$true_shared_factors),baseline$true_spfac)
})
