## -----------------------------------------------------------------------------
library(MSGPCA)
s <- simulation_gg_blocks(nside=6,Jsp=8,Jmix=12,Jns=10,seed=101)
f <- MSGPCA(lapply(s$Y_list,log1p),s$locs_list,n_shared=3,n_local=1,
             lambda_shared=c(.1,.3,.6),lambda_local=c(.1,.3),
             tuning='pareto',seed=101,verbose=FALSE)
f$optimal_lambda_shared
f$optimal_lambda_local
f$pareto_results
plot_pareto(f)

## ----eval=FALSE---------------------------------------------------------------
# # Full empirical grid (kept out of routine vignette checks).
# f <- tune_MSGPCA(lapply(s$Y_list,log1p),s$locs_list,
#                  lambda_shared=seq(.1,.6,.1),lambda_local=seq(.1,.6,.1),
#                  n_shared=3,n_local=1,lambda_constraint='shared_ge_local',
#                  seed=101,verbose=FALSE)

