## -----------------------------------------------------------------------------
library(MSGPCA)
sim <- simulation_gg_blocks(nside=12,seed=101)
names(sim)
dim(sim$Y_list[[1]])
sim$gene_groups
sim$perturbations[c('affected_slices','intensity','mechanism')]

## -----------------------------------------------------------------------------
slab <- simulation_gg_blocks(nside=6,scenario='slab',seed=101)
none <- simulation_gg_blocks(nside=6,scenario='none',seed=101)
cohort <- simulation_gg_blocks(nside=6,n_slices=5,
                              scenario='centroid_shift',seed=101)
length(cohort$true_cluster_labels)

